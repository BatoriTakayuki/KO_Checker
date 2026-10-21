#!/usr/bin/env python3
"""RNA evidence for a reported KO. Coordinates in CSVs are genomic 1-based."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import urllib.request
from transfers import download_file

try:
    import pysam
except ImportError:
    pysam = None


REFERENCE = {
    "mm10": ("Gencode_mouse/release_M25", "GRCm38.primary_assembly.genome.fa.gz", "gencode.vM25.annotation.gtf.gz"),
    "mm39": ("Gencode_mouse/release_M39", "GRCm39.primary_assembly.genome.fa.gz", "gencode.vM39.annotation.gtf.gz"),
    "hg38": ("Gencode_human/release_47", "GRCh38.primary_assembly.genome.fa.gz", "gencode.v47.annotation.gtf.gz"),
}
BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/"
CODONS = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
TRANS = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.upper().translate(TRANS)[::-1]


def translate(dna):
    protein = ""
    stop_offset = None
    for i in range(0, len(dna) - 2, 3):
        aa = CODONS.get(dna[i:i+3].upper(), "X")
        if aa == "*":
            stop_offset = i
            break
        protein += aa
    return protein, stop_offset


def attrs_parse(s):
    return dict(re.findall(r'(\w+) "([^"]*)"', s))


@dataclass
class Transcript:
    tid: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    biotype: str = ""
    exons: list = field(default_factory=list)
    cds: list = field(default_factory=list)

    def ordered_exons(self):
        return sorted(set(self.exons), reverse=(self.strand == "-"))

    def ordered_cds(self):
        return sorted(set(self.cds), reverse=(self.strand == "-"))


def load_gene(gtf, gene):
    opener = gzip.open if str(gtf).endswith(".gz") else open
    tx = {}
    folded_tx = {}
    with opener(gtf, "rt") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            bits = line.rstrip("\n").split("\t")
            if len(bits) != 9 or bits[2] not in ("exon", "CDS"):
                continue
            attrs = attrs_parse(bits[8])
            exact = gene in (attrs.get("gene_name"), attrs.get("gene_id"))
            if not exact and gene.casefold() != attrs.get("gene_name", "").casefold():
                continue
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            gene_id = attrs.get("gene_id", "")
            matches = tx if exact else folded_tx
            t = matches.setdefault(tid, Transcript(tid, gene_id, attrs.get("gene_name", ""), bits[0], bits[6], attrs.get("gene_type", attrs.get("gene_biotype", ""))))
            (t.exons if bits[2] == "exon" else t.cds).append((int(bits[3]) - 1, int(bits[4])))
    tx = tx or folded_tx
    if not tx:
        raise ValueError(f"Gene {gene!r} not found in GTF")
    genes = {t.gene_id for t in tx.values()}
    if len(genes) != 1:
        raise ValueError(f"Gene name {gene!r} is ambiguous; supply a gene_id instead ({sorted(genes)})")
    return {k: v for k, v in tx.items() if v.exons}


def download(url, dest, progress=None, check_pause=None):
    download_file(url, dest, progress, check_pause)


def prepare_reference(assembly, cache, fasta=None, gtf=None, progress=None, check_pause=None):
    progress = progress or (lambda **kw: None)
    check_pause = check_pause or (lambda: None)
    check_pause()
    if fasta and gtf:
        fa, annotation = Path(fasta), Path(gtf)
    elif fasta or gtf:
        raise ValueError("Supply both --fasta and --gtf, or neither")
    else:
        directory, fa_name, gtf_name = REFERENCE[assembly]
        cache = Path(cache) / assembly
        cache.mkdir(parents=True, exist_ok=True)
        compressed = cache / fa_name
        annotation = cache / gtf_name
        fa = cache / fa_name[:-3]
        for path in (compressed, annotation):
            if not path.exists():
                progress(stage='reference_download', message=f'参照を取得中: {path.name}', transfer=None)
                download(BASE_URL + directory + "/" + path.name, path,
                         lambda **kw: progress(transfer=kw), check_pause)
            else:
                progress(stage='reference_download', message=f'参照を再利用: {path.name}', transfer=None)
        if not fa.exists():
            progress(stage='reference_extract', message='参照ゲノムを展開中', transfer=None)
            tmp = fa.with_name(fa.name + ".partial")
            try:
                with gzip.open(compressed, "rb") as src, open(tmp, "wb") as dst:
                    import time
                    notified = 0
                    while True:
                        check_pause()
                        chunk = src.read(1024 * 1024)
                        if not chunk:
                            break
                        dst.write(chunk)
                        if time.monotonic() - notified >= 1:
                            progress(activity=f'展開済み {dst.tell() / 1000000:.1f} MB')
                            notified = time.monotonic()
                os.replace(tmp, fa)
            finally:
                tmp.unlink(missing_ok=True)
    for path in (fa, annotation):
        if not path.exists():
            raise FileNotFoundError(str(path))
    if not Path(str(fa) + ".fai").exists():
        check_pause()
        progress(stage='reference_index', message='FASTAの索引を作成中', transfer=None)
        pysam.faidx(str(fa))
    check_pause()
    return fa, annotation


def read_samples(args):
    records = []
    if args.samples:
        with open(args.samples, newline="") as f:
            records = list(csv.DictReader(f))
    else:
        # CLI: --wt SRX1,SRX2 --ko SRX3; metadata labels may also be KD/OE.
        for label in ("WT", "KO", "KD", "OE"):
            value = getattr(args, label.lower())
            records += [{"srx_id": s.strip(), "metadata_label": label} for s in value.split(",") if s.strip()]
    if not records:
        raise ValueError("Provide --samples or at least one of --wt/--ko/--kd/--oe")
    seen = set()
    for record in records:
        sid = record.get("srx_id", "").strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sid) or sid in seen:
            raise ValueError(f"Invalid or duplicate srx_id: {sid!r}")
        seen.add(sid)
        label = record.get("metadata_label", "").strip().upper()
        if label not in ("WT", "KO", "KD", "OE"):
            raise ValueError(f"Unknown metadata_label for {sid}: {label!r}")
        record["srx_id"], record["metadata_label"] = sid, label
        path = record.get("bam_path", "").strip()
        if not path:
            path = args.bam_pattern.format(srx=sid)
        record["bam_path"] = path
        if not Path(path).is_file():
            raise FileNotFoundError(f"BAM for {sid}: {path}")
    return records


def resolve_reference_name(chrom, reference):
    variants = [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]
    for name in variants:
        if name in reference.references:
            return name
    raise ValueError(f"Contig {chrom} is not present in reference")


def resolve_contig(chrom, bam, fa):
    bam_name = resolve_reference_name(chrom, bam)
    fasta_name = resolve_reference_name(chrom, fa)
    if bam.get_reference_length(bam_name) != fa.get_reference_length(fasta_name):
        raise ValueError(f"BAM/FASTA contig length differs for {bam_name}/{fasta_name}; check reference assembly")
    return bam_name


def read_ok(read, min_mapq):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
        return False
    if read.has_tag("NH"):
        return read.get_tag("NH") == 1 and read.mapping_quality >= min_mapq
    return read.mapping_quality >= max(min_mapq, 20)


def exonic_overlap(read, exons):
    # Only actual aligned blocks; introns in a spliced read are excluded.
    return any(a < e and b > s for a, b in read.get_blocks() for s, e in exons)


def gene_exons(tx):
    intervals = sorted({x for t in tx.values() for x in t.exons})
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
        else:
            merged.append((s, e))
    return merged


def junctions_from_read(read):
    pos = read.reference_start
    for op, length in read.cigartuples or []:
        if op == 3:
            yield (pos, pos + length)
        if op in (0, 2, 3, 7, 8):
            pos += length


def count_library(bam, min_mapq):
    # Scan once. read1 for paired ends, all unpaired reads for single ends.
    n = 0
    for r in bam.fetch(until_eof=True):
        if read_ok(r, min_mapq) and (not r.is_paired or r.is_read1):
            n += 1
    return n


def collect(bam, chrom, start, end, exons, min_mapq):
    fragments = set()
    junction_fragments = defaultdict(set)
    for read in bam.fetch(chrom, start, end):
        if not read_ok(read, min_mapq):
            continue
        for junction in junctions_from_read(read):
            junction_fragments[junction].add(read.query_name)
        if exonic_overlap(read, exons):
            fragments.add(read.query_name)
    return len(fragments), {k: len(v) for k, v in junction_fragments.items()}


def junction_observation(junctions, junction):
    """Local junction support, using the better-covered shared boundary.

    Counts are junction-support fragments, not an isoform abundance estimate.
    Taking the larger boundary total avoids double-counting two normal junctions.
    """
    left, right = junction
    depth = max(sum(n for (a, b), n in junctions.items() if a == left),
                sum(n for (a, b), n in junctions.items() if b == right))
    return depth, junctions.get(junction, 0)


def call_variants(bam, fa, chrom, start, end, exons, min_mapq, min_bq, min_alt, min_af):
    calls = []
    fasta_chrom = resolve_reference_name(chrom, fa)
    for col in bam.pileup(chrom, start, end, truncate=True, stepper="nofilter", min_base_quality=0, max_depth=1000000, ignore_overlaps=True):
        pos = col.reference_pos
        if not any(s <= pos < e for s, e in exons):
            continue
        reference_base = fa.fetch(fasta_chrom, pos, pos + 1).upper()
        if reference_base not in "ACGT":
            continue
        alleles = Counter()
        strand_support = defaultdict(lambda: [0, 0])
        for p in col.pileups:
            r, q = p.alignment, p.query_position
            if not read_ok(r, min_mapq) or p.is_refskip or p.is_del or q is None:
                continue
            if r.query_qualities is None or r.query_qualities[q] < min_bq:
                continue
            base = r.query_sequence[q].upper()
            if base not in "ACGT":
                continue
            key = (reference_base, base)
            if p.indel > 0:
                ins = r.query_sequence[q+1:q+1+p.indel].upper()
                if len(ins) != p.indel or set(ins) - set("ACGT"):
                    continue
                key = (reference_base, reference_base + ins)
            elif p.indel < 0:
                deleted = fa.fetch(fasta_chrom, pos+1, pos+1-p.indel).upper()
                key = (reference_base + deleted, reference_base)
            alleles[key] += 1
            strand_support[key][int(r.is_reverse)] += 1
        depth = sum(alleles.values())
        if not depth:
            continue
        for (ref, alt), count in alleles.items():
            if ref == alt or count < min_alt or count / depth < min_af:
                continue
            calls.append({"chrom": chrom, "pos": pos + 1, "ref": ref, "alt": alt,
                          "kind": "SNV" if len(ref) == len(alt) == 1 else "INDEL",
                          "alt_reads": count, "depth": depth, "allele_fraction": round(count / depth, 5),
                          "forward_reads": strand_support[(ref, alt)][0], "reverse_reads": strand_support[(ref, alt)][1]})
    return calls


def observe_event(bam, chrom, event, min_mapq, min_bq):
    """Recheck a candidate in every sample, including samples below calling thresholds."""
    pos = int(event["pos_1based"]) - 1
    depth = alt_count = 0
    for col in bam.pileup(chrom, pos, pos + 1, truncate=True, stepper="nofilter", min_base_quality=0, max_depth=1000000, ignore_overlaps=True):
        if col.reference_pos != pos:
            continue
        for p in col.pileups:
            r, q = p.alignment, p.query_position
            if not read_ok(r, min_mapq) or p.is_refskip or p.is_del or q is None:
                continue
            if r.query_qualities is None or r.query_qualities[q] < min_bq:
                continue
            base = r.query_sequence[q].upper()
            if base not in "ACGT":
                continue
            # A usable anchor is the denominator, irrespective of its allele.
            depth += 1
            ref, alt = event["ref"], event["alt"]
            if len(ref) == len(alt) == 1:
                is_alt = base == alt
            elif len(alt) > len(ref):
                insertion = r.query_sequence[q+1:q+1+p.indel].upper() if p.indel > 0 else ""
                is_alt = base == ref[0] and insertion == alt[len(ref):]
            else:
                is_alt = base == ref[0] and p.indel == len(alt) - len(ref)
            alt_count += int(is_alt)
    return depth, alt_count


def template(t, fa, chrom):
    exons = t.ordered_exons()
    seq, g2t = "", {}
    for start, end in exons:
        part = fa.fetch(chrom, start, end).upper()
        if t.strand == "-":
            part = rc(part)
            positions = range(end - 1, start - 1, -1)
        else:
            positions = range(start, end)
        offset = len(seq)
        seq += part
        for idx, gpos in enumerate(positions):
            g2t[gpos] = offset + idx
    cds_idxs = [g2t[pos] for s, e in t.ordered_cds() for pos in range(s, e) if pos in g2t]
    cds_start = min(cds_idxs) if cds_idxs else None
    cds_end = max(cds_idxs) + 1 if cds_idxs else None
    boundaries = [sum(e - s for s, e in exons[:i]) for i in range(1, len(exons))]
    return seq, g2t, cds_start, cds_end, boundaries


def altered_transcript(seq, g2t, t, event):
    pos = event["pos"] - 1
    ref, alt = event["ref"], event["alt"]
    if len(ref) == len(alt) == 1:
        if pos not in g2t:
            return None
        idx = g2t[pos]
        original, change = (rc(ref), rc(alt)) if t.strand == "-" else (ref, alt)
        if seq[idx] != original:
            return None
        return seq[:idx] + change + seq[idx+1:], idx, 1, 1
    if len(alt) > len(ref) and alt.startswith(ref):
        if pos not in g2t:
            return None
        idx = g2t[pos]
        inserted = alt[len(ref):] if t.strand == "+" else rc(alt[len(ref):])
        insert_at = idx + 1 if t.strand == "+" else idx
        return seq[:insert_at] + inserted + seq[insert_at:], insert_at, 0, len(inserted)
    if len(ref) > len(alt) and ref.startswith(alt):
        deleted = list(range(pos + len(alt), pos + len(ref)))
        if not all(p in g2t for p in deleted):
            return None
        offsets = sorted(g2t[p] for p in deleted)
        if offsets != list(range(offsets[0], offsets[-1] + 1)):
            return None
        idx = offsets[0]
        expected = rc(ref[len(alt):]) if t.strand == "-" else ref[len(alt):]
        if seq[idx:idx+len(expected)] != expected:
            return None
        return seq[:idx] + seq[idx+len(deleted):], idx, len(deleted), 0
    return None


def protein_effect(t, event, fa, chrom):
    seq, g2t, cds_start, cds_end, boundaries = template(t, fa, chrom)
    out = {"transcript_id": t.tid, "gene_biotype": t.biotype,
           "effect": "noncoding_or_unresolved", "reference_protein": "", "predicted_protein": "",
           "first_changed_aa": "", "stop_codon_cdna_1based": "", "stop_to_last_junction_nt": "", "nmd_prediction": "unknown",
           "protein_prediction_status": "not_available"}
    if cds_start is None or seq[cds_start:cds_start+3] != "ATG":
        return out
    ref_protein, ref_stop = translate(seq[cds_start:])
    out["reference_protein"] = ref_protein
    mutation = altered_transcript(seq, g2t, t, event)
    if mutation is None:
        out["effect"] = "outside_exon_or_complex_indel"
        return out
    new_seq, edit_idx, removed, added = mutation
    if edit_idx < cds_start or edit_idx >= cds_end:
        out["effect"] = "noncoding_exon_or_utr"
        return out
    if edit_idx < cds_start + 3:
        out["effect"] = "start_codon_affected"
        return out
    shifted_cds = new_seq[cds_start:]
    new_protein, stop_offset = translate(shifted_cds)
    out["predicted_protein"] = new_protein
    out["protein_prediction_status"] = "single_event_model; other distant edits not phased"
    delta = added - removed
    if delta % 3:
        out["effect"] = "frameshift"
    elif ref_protein == new_protein:
        out["effect"] = "synonymous_or_same_predicted_protein"
    elif stop_offset is not None and ref_stop is not None and stop_offset < ref_stop + delta:
        out["effect"] = "stop_gained"
    elif delta < 0:
        out["effect"] = "inframe_deletion"
    elif delta > 0:
        out["effect"] = "inframe_insertion"
    else:
        out["effect"] = "missense_or_inframe_change"
    first = next((i + 1 for i in range(min(len(ref_protein), len(new_protein))) if ref_protein[i] != new_protein[i]), None)
    if first is None and len(ref_protein) != len(new_protein):
        first = min(len(ref_protein), len(new_protein)) + 1
    out["first_changed_aa"] = first or ""
    if stop_offset is not None:
        stop_pos = cds_start + stop_offset
        out["stop_codon_cdna_1based"] = stop_pos + 1
        if boundaries:
            last_boundary = boundaries[-1] + (delta if edit_idx < boundaries[-1] else 0)
            distance = last_boundary - (stop_pos + 3)
            out["stop_to_last_junction_nt"] = distance
            out["nmd_prediction"] = "possible" if distance > 55 else "unlikely_by_50_55nt_rule" if distance <= 50 else "borderline"
        else:
            out["nmd_prediction"] = "unlikely_single_exon"
    return out


def expected_junctions(t):
    exons = t.ordered_exons()
    return {(min(a[1], b[1]), max(a[0], b[0])) for a, b in zip(exons, exons[1:])}


def skipped_exons(t, junction):
    exons = t.ordered_exons()
    for i, a in enumerate(exons):
        for j in range(i+2, len(exons)):
            b = exons[j]
            candidate = (min(a[1], b[1]), max(a[0], b[0]))
            if candidate == junction:
                return exons[i+1:j]
    return []


def junction_effect(t, junction, fa, chrom):
    skipped = skipped_exons(t, junction)
    if not skipped:
        return None
    seq, g2t, cds_start, cds_end, _ = template(t, fa, chrom)
    if cds_start is None or seq[cds_start:cds_start+3] != "ATG":
        return {"effect": "exon_skip_noncoding_or_unresolved", "reference_protein": "", "predicted_protein": "", "stop_codon_cdna_1based": "", "stop_to_last_junction_nt": "", "nmd_prediction": "unknown", "first_changed_aa": "", "protein_prediction_status": "not_available"}
    positions = sorted(g2t[p] for s, e in skipped for p in range(s, e) if p in g2t)
    if not positions or positions != list(range(positions[0], positions[-1]+1)):
        return None
    first, last = positions[0], positions[-1]+1
    if first <= cds_start < last:
        return {"effect": "start_region_skipped", "reference_protein": translate(seq[cds_start:])[0], "predicted_protein": "", "stop_codon_cdna_1based": "", "stop_to_last_junction_nt": "", "nmd_prediction": "unknown", "first_changed_aa": "", "protein_prediction_status": "not_available"}
    new_seq = seq[:first] + seq[last:]
    new_cds_start = cds_start - (last-first if last <= cds_start else 0)
    reference_protein, _ = translate(seq[cds_start:])
    predicted, stop = translate(new_seq[new_cds_start:])
    delta = last-first
    surviving_exons = [ex for ex in t.ordered_exons() if ex not in skipped]
    remaining_boundaries = [sum(e-s for s,e in surviving_exons[:i]) for i in range(1,len(surviving_exons))]
    last_boundary = max(remaining_boundaries) if remaining_boundaries else None
    stop_pos = new_cds_start + stop if stop is not None else None
    distance = last_boundary-(stop_pos+3) if last_boundary is not None and stop_pos is not None else None
    first_aa = next((i+1 for i in range(min(len(reference_protein),len(predicted))) if reference_protein[i]!=predicted[i]), min(len(reference_protein),len(predicted))+1)
    effect = "exon_skip_outside_cds" if reference_protein == predicted else "exon_skip_frameshift" if delta % 3 else "exon_skip_inframe"
    return {"effect": effect,
            "reference_protein":reference_protein,"predicted_protein":predicted,"first_changed_aa":first_aa,
            "stop_codon_cdna_1based":stop_pos+1 if stop_pos is not None else "",
            "stop_to_last_junction_nt": distance if distance is not None else "",
            "nmd_prediction":("possible" if distance>55 else "unlikely_by_50_55nt_rule" if distance<=50 else "borderline") if distance is not None else "unknown",
            "protein_prediction_status":"single_skip_model; other distant edits not phased"}


def write_csv(path, rows, columns):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    if pysam is None:
        raise RuntimeError("pysam is required (use the supplied Docker image)")
    if args.min_alt < 1 or not 0 <= args.min_af <= 1 or args.min_junction < 1:
        raise ValueError("Thresholds must be positive, min-af between 0 and 1")
    samples = read_samples(args)
    fasta_path, gtf = prepare_reference(args.assembly, args.cache, args.fasta, args.gtf)
    tx = load_gene(gtf, args.gene)
    chromosomes = {t.chrom for t in tx.values()}
    if len(chromosomes) != 1:
        raise ValueError("Gene transcripts have inconsistent chromosomes")
    annotated_chrom = chromosomes.pop()
    exons = gene_exons(tx)
    start, end = exons[0][0], exons[-1][1]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    expected = set().union(*(expected_junctions(t) for t in tx.values()))
    fa = pysam.FastaFile(str(fasta_path))
    fasta_chrom = resolve_reference_name(annotated_chrom, fa)
    sample_rows, event_rows, effect_rows, junction_rows = [], [], [], []
    sample_junctions = {}
    for sample in samples:
        sid, label = sample["srx_id"], sample["metadata_label"]
        with pysam.AlignmentFile(sample["bam_path"], "rb") as bam:
            if not bam.has_index():
                raise ValueError(f"BAM index missing for {sid}; create .bam.bai before analysis")
            chrom = resolve_contig(annotated_chrom, bam, fa)
            library_fragments = count_library(bam, args.min_mapq)
            gene_fragments, junctions = collect(bam, chrom, start, end, exons, args.min_mapq)
            calls = call_variants(bam, fa, chrom, start, end, exons, args.min_mapq, args.min_bq, args.min_alt, args.min_af)
        cpm = round(gene_fragments / library_fragments * 1000000, 5) if library_fragments else ""
        sample_rows.append({"srx_id":sid,"metadata_label":label,"bam_path":sample["bam_path"],
                            "gene":args.gene,"assembly":args.assembly,"reference_gtf":str(gtf),
                            "unique_primary_library_fragments":library_fragments,"unique_exonic_gene_fragments":gene_fragments,
                            "gene_unique_fragment_cpm":cpm,"wt_baseline_cpm":"","expression_relative_to_wt":"",
                            "predicted_class":"undetermined","classification_reason":"","coverage_warning":"low_exonic_fragments" if gene_fragments < args.min_gene_fragments else "",
                            "normal_mrna_fraction":"","abnormal_mrna_fraction":"","mrna_fraction_reason":"distant variants and isoforms not phased; local allele fractions are in events.csv"})
        for call in calls:
            eid = hashlib.sha1(f"{args.assembly}:{call['chrom']}:{call['pos']}:{call['ref']}:{call['alt']}".encode()).hexdigest()[:16]
            row = {"event_id":eid,"srx_id":sid,"metadata_label":label,"gene":args.gene,"event_type":call["kind"],
                   "chrom":chrom,"pos_1based":call["pos"],"ref":call["ref"],"alt":call["alt"],
                   "alt_reads":call["alt_reads"],"depth_at_anchor":call["depth"],"local_alt_fraction":call["allele_fraction"],
                   "forward_reads":call["forward_reads"],"reverse_reads":call["reverse_reads"],"wt_label_max_fraction":"","comparison_status":""}
            event_rows.append(row)
            for t in tx.values():
                if not any(s <= call["pos"]-1 < e for s,e in t.exons):
                    continue
                effect_rows.append({"event_id":eid,"srx_id":sid,"metadata_label":label,"gene":args.gene,
                                    **protein_effect(t, call, fa, fasta_chrom)})
        sample_junctions[sid] = junctions
        for (jstart, jend), n in junctions.items():
            if n < args.min_junction:
                continue
            jid = hashlib.sha1(f"{args.assembly}:{chrom}:{jstart}:{jend}:splice".encode()).hexdigest()[:16]
            skip_tx = [t for t in tx.values() if skipped_exons(t,(jstart,jend))]
            junction_rows.append({"event_id":jid,"srx_id":sid,"metadata_label":label,"gene":args.gene,"chrom":chrom,
                                  "intron_start_1based":jstart+1,"intron_end_1based":jend,"support_fragments":n,
                                  "annotation_status":"annotated" if (jstart,jend) in expected else "unannotated",
                                  "skipped_transcripts":";".join(t.tid for t in skip_tx),"shared_with_wt_label":""})
            for t in skip_tx:
                consequence = junction_effect(t,(jstart,jend),fa,fasta_chrom)
                if consequence:
                    effect_rows.append({"event_id":jid,"srx_id":sid,"metadata_label":label,"gene":args.gene,
                                        "transcript_id":t.tid,"gene_biotype":t.biotype,**consequence})
    wt_samples = {r["srx_id"] for r in sample_rows if r["metadata_label"] == "WT"}
    distinct_events = {r["event_id"]: r for r in event_rows}
    observations = {}
    for s in sample_rows:
        with pysam.AlignmentFile(s["bam_path"], "rb") as bam:
            chrom = resolve_contig(annotated_chrom, bam, fa)
            for eid, candidate in distinct_events.items():
                observations[(s["srx_id"], eid)] = observe_event(bam, chrom, candidate, args.min_mapq, args.min_bq)
    by_sample_events = defaultdict(list)
    for r in event_rows:
        candidates = [sid for sid in wt_samples if sid != r["srx_id"]]
        if not candidates:
            candidates = [s["srx_id"] for s in sample_rows if s["srx_id"] != r["srx_id"]]
        covered = [(alt / depth, sid) for sid in candidates
                   for depth, alt in [observations[(sid, r["event_id"])]] if depth >= args.min_event_depth]
        max_other = max((fraction for fraction, _ in covered), default=None)
        r["wt_label_max_fraction"] = round(max_other, 5) if max_other is not None else ""
        r["comparison_status"] = ("enriched_vs_comparators" if max_other is not None and float(r["local_alt_fraction"]) - max_other >= args.min_delta
                                  else "shared_or_not_enriched" if max_other is not None else "no_covered_comparator")
        by_sample_events[r["srx_id"]].append(r)
    skip_event_map = {}
    for r in junction_rows:
        junction = (int(r["intron_start_1based"])-1, int(r["intron_end_1based"]))
        skip_event_map[r["event_id"]] = junction
        comparators = [sid for sid in wt_samples if sid != r["srx_id"]]
        if not comparators:
            comparators = [sid for sid in sample_junctions if sid != r["srx_id"]]
        r["shared_with_wt_label"] = ("yes" if any(sample_junctions[sid].get(junction, 0) >= args.min_junction for sid in comparators)
                                     else "no" if comparators else "unknown")
        depth, support = junction_observation(sample_junctions[r["srx_id"]], junction)
        covered = [n / d for sid in comparators
                   for d, n in [junction_observation(sample_junctions[sid], junction)]
                   if d >= args.min_event_depth]
        fraction = support / depth if depth else 0
        maximum = max(covered) if covered else None
        r["local_junction_depth"] = depth
        r["local_skip_fraction"] = round(fraction, 5)
        r["comparator_max_skip_fraction"] = round(maximum, 5) if maximum is not None else ""
        r["comparison_status"] = ("enriched_vs_comparators" if depth >= args.min_event_depth
                                  and maximum is not None and fraction - maximum >= args.min_delta
                                  else "shared_or_not_enriched" if maximum is not None else "no_covered_comparator")
    damaging = {"frameshift","stop_gained","exon_skip_frameshift","start_codon_affected","start_region_skipped"}
    ko_evidence = {}
    ko_skip_evidence = {}
    for s in sample_rows:
        sid = s["srx_id"]
        enriched = {r["event_id"] for r in by_sample_events[sid] if r["comparison_status"] == "enriched_vs_comparators" and float(r["local_alt_fraction"]) >= args.min_ko_af}
        detrimental = {r["event_id"] for r in effect_rows if r["srx_id"] == sid and r["effect"] in damaging}
        coding_skips = {r['event_id'] for r in effect_rows if r['srx_id'] == sid
                        and r['effect'] in ('exon_skip_frameshift', 'exon_skip_inframe', 'start_region_skipped')}
        novel_skip = {r["event_id"] for r in junction_rows if r["srx_id"] == sid and r["event_id"] in coding_skips and r["comparison_status"] == "enriched_vs_comparators" and float(r["local_skip_fraction"]) >= args.min_ko_af and r["annotation_status"] == "unannotated"}
        ko_evidence[sid] = enriched & detrimental
        ko_skip_evidence[sid] = novel_skip
        if s["coverage_warning"]:
            s["classification_reason"] = "Insufficient uniquely mapped exonic fragments to interpret lack of evidence"
        elif (enriched & detrimental) or novel_skip:
            s["predicted_class"] = "KO_like"
            s["classification_reason"] = "Potential damaging event enriched relative to WT-labelled samples; RNA evidence only"
        else:
            s["classification_reason"] = "No KO-specific disruptive event demonstrated; cannot establish wild-type genotype from RNA"
    # A WT-like call requires enough reference observations at a KO-like variant.
    ko_like_events = set().union(*(ko_evidence[s["srx_id"]] for s in sample_rows if s["predicted_class"] == "KO_like"))
    ko_like_skips = set().union(*(ko_skip_evidence[s["srx_id"]] for s in sample_rows if s["predicted_class"] == "KO_like"))
    for s in sample_rows:
        sid = s["srx_id"]
        if s["predicted_class"] != "undetermined" or s["coverage_warning"]:
            continue
        checked = [(depth, alt) for eid in ko_like_events for depth, alt in [observations[(s["srx_id"], eid)]]
                   if depth >= args.min_event_depth]
        normal_splicing = False
        for eid in ko_like_skips:
            junction = skip_event_map[eid]
            splice_depth, skip_support = junction_observation(sample_junctions[sid], junction)
            if splice_depth < args.min_event_depth or skip_support / splice_depth >= args.min_af:
                continue
            for t in tx.values():
                skipped = skipped_exons(t, junction)
                skip_boundaries = {v for start, end in skipped for v in (start, end)}
                if skipped and any(sample_junctions[sid].get(j, 0) >= args.min_junction
                                   for j in expected_junctions(t) if j[0] in skip_boundaries or j[1] in skip_boundaries):
                    normal_splicing = True
                    break
        if (checked or normal_splicing) and all(alt/depth < args.min_af for depth, alt in checked):
            s["predicted_class"] = "WT_like"
            s["classification_reason"] = "Reference bases at covered KO-like variant site(s) and/or annotated junction support; other lesions remain possible"
    # Expression comparison is independent of sequence-based classification.
    candidate_wt = [s for s in sample_rows if s["metadata_label"] == "WT" and s["gene_unique_fragment_cpm"] != ""]
    baseline = sum(s["gene_unique_fragment_cpm"] for s in candidate_wt)/len(candidate_wt) if candidate_wt else None
    for s in sample_rows:
        if baseline and baseline > 0:
            s["wt_baseline_cpm"] = round(baseline,5)
            s["expression_relative_to_wt"] = round(s["gene_unique_fragment_cpm"]/baseline,5) if s["gene_unique_fragment_cpm"] != "" else ""
    sample_cols = ["srx_id","metadata_label","bam_path","gene","assembly","reference_gtf","unique_primary_library_fragments","unique_exonic_gene_fragments","gene_unique_fragment_cpm","wt_baseline_cpm","expression_relative_to_wt","predicted_class","classification_reason","coverage_warning","normal_mrna_fraction","abnormal_mrna_fraction","mrna_fraction_reason"]
    event_cols = ["event_id","srx_id","metadata_label","gene","event_type","chrom","pos_1based","ref","alt","alt_reads","depth_at_anchor","local_alt_fraction","forward_reads","reverse_reads","wt_label_max_fraction","comparison_status"]
    junction_cols = ["event_id","srx_id","metadata_label","gene","chrom","intron_start_1based","intron_end_1based","support_fragments","annotation_status","skipped_transcripts","shared_with_wt_label","local_junction_depth","local_skip_fraction","comparator_max_skip_fraction","comparison_status"]
    effect_cols = ["event_id","srx_id","metadata_label","gene","transcript_id","gene_biotype","effect","reference_protein","predicted_protein","first_changed_aa","stop_codon_cdna_1based","stop_to_last_junction_nt","nmd_prediction","protein_prediction_status"]
    for name, rows, cols in (("samples.csv",sample_rows,sample_cols),("events.csv",event_rows,event_cols),("junctions.csv",junction_rows,junction_cols),("consequences.csv",effect_rows,effect_cols)):
        write_csv(out_dir/name,rows,cols)
    manifest = {"version":"0.1.2","assembly":args.assembly,"gene":args.gene,"fasta":str(fasta_path),"gtf":str(gtf),
                "expression_baseline":"all_WT_labelled_samples_with_valid_library",
                "confirmed_replicates":getattr(args, 'confirm_replicates', False),
                "thresholds":{"min_alt":args.min_alt,"min_af":args.min_af,"min_junction":args.min_junction,"min_delta":args.min_delta,"min_gene_fragments":args.min_gene_fragments,"min_event_depth":args.min_event_depth,"min_ko_af":args.min_ko_af,"min_mapq":args.min_mapq,"min_bq":args.min_bq},
                "notes":["RNA evidence is not a DNA KO validation","CPM counts uniquely mapped primary fragments, not TPM","Variant fractions are local and not whole-mRNA fractions","Predicted proteins model one event at a time"]}
    (out_dir/"run.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    fa.close()
    if not getattr(args, 'skip_region_analysis', False):
        from pipeline import write_exon_coverage, write_observations
        from region_loss import write_region_loss
        write_observations(samples, out_dir, fasta_path, gtf, args.gene)
        write_exon_coverage(samples, out_dir, fasta_path, gtf, args.gene)
        write_region_loss(out_dir, fasta_path, gtf, args.gene)
    from reporting import write_conclusion
    write_conclusion(out_dir, getattr(args, 'confirm_replicates', False))
    print(f"Finished: {out_dir}")


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene",required=True,help="Gene symbol or stable gene ID")
    p.add_argument("--assembly",choices=REFERENCE,default="mm10")
    p.add_argument("--samples",help="CSV columns: srx_id,metadata_label,bam_path (path optional)")
    p.add_argument('--skip-region-analysis', action='store_true', help=argparse.SUPPRESS)
    p.add_argument('--confirm-replicates', action='store_true',
                   help='Confirm independent biological replicates and exchangeable, comparable conditions for an exploratory CPM permutation test')
    for label in ("wt","ko","kd","oe"):
        p.add_argument("--"+label,default="",help="Comma-separated SRX IDs")
    p.add_argument("--bam-pattern",default="/data/Analysis/{srx}/{srx}.bam")
    p.add_argument("--fasta",help="Override downloaded reference FASTA; requires --gtf")
    p.add_argument("--gtf",help="Override annotation GTF; requires --fasta")
    p.add_argument("--cache",default="/references")
    p.add_argument("--out",default="/results")
    p.add_argument("--min-mapq",type=int,default=0,help="NH=1 required when available; otherwise MAPQ>=20")
    p.add_argument("--min-bq",type=int,default=20)
    p.add_argument("--min-alt",type=int,default=3)
    p.add_argument("--min-af",type=float,default=0.05)
    p.add_argument("--min-junction",type=int,default=3)
    p.add_argument("--min-delta",type=float,default=0.2)
    p.add_argument("--min-ko-af",type=float,default=0.4)
    p.add_argument("--min-event-depth",type=int,default=10)
    p.add_argument("--min-gene-fragments",type=int,default=20)
    return p


if __name__ == "__main__":
    try:
        run(parser().parse_args())
    except (ValueError,FileNotFoundError,RuntimeError,OSError) as e:
        print(f"ERROR: {e}",file=sys.stderr)
        sys.exit(2)
