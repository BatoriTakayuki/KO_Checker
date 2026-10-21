"""Infer RNA exon-region loss and explicitly conditional transcript/protein models.

Coverage does not establish a DNA deletion or an unobserved junction. Compare each
exon to library-normalized WT, requiring a well-covered retained region as an
internal control instead of requiring the missing exons themselves to be covered.
"""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
from statistics import mean

import rnako


REASON = 'KO-like RNA region loss in a major coding transcript; DNA breakpoints and protein function unconfirmed'
REFERENCE_REASON = 'Reference-like RNA coverage across candidate lost coding regions; DNA genotype unconfirmed'


def apply_region_classes(samples, exons, models, min_depth=10):
    """Assess local region retention independently of the supplied sample label."""
    supported = {m['srx_id'] for m in models if m['ko_support'] == 'yes'}
    regions = set()
    for model in models:
        if model['ko_support'] == 'yes':
            for region in model['lost_regions'].split(';'):
                chrom, interval = region.rsplit(':', 1)
                start, end = interval.split('-')
                regions.add((chrom, int(start), int(end)))
    coverage = defaultdict(dict)
    for row in exons:
        if not row.get('chrom'):
            continue
        coverage[row['srx_id']][(row['chrom'], int(row['exon_start_1based']), int(row['exon_end_1based']))] = row
    for sample in samples:
        if sample.get('classification_reason') in (REASON, REFERENCE_REASON):
            sample.update(predicted_class='undetermined', classification_reason='No KO-specific disruptive event demonstrated; cannot establish wild-type genotype from RNA')
        if sample['srx_id'] in supported and sample['predicted_class'] != 'KO_like':
            sample.update(predicted_class='KO_like', classification_reason=REASON)
        elif (regions and sample['predicted_class'] == 'undetermined'
              and not sample.get('coverage_warning')
              and float(sample.get('unique_primary_library_fragments') or 0) > 0):
            rows = coverage[sample['srx_id']]
            if all(region in rows and float(rows[region]['mean_depth']) >= min_depth
                   and float(rows[region]['sufficient_depth_fraction']) >= .8 for region in regions):
                sample.update(predicted_class='WT_like', classification_reason=REFERENCE_REASON)
COLUMNS = ['event_id', 'srx_id', 'metadata_label', 'transcript_id', 'chrom', 'strand',
           'lost_exons', 'lost_regions', 'wt_reference_samples', 'maximum_normalized_depth_ratio',
           'retained_control_regions', 'evidence', 'junction_support', 'ko_support', 'effect',
           'deleted_coding_nt', 'coding_loss_fraction', 'reference_protein', 'predicted_protein',
           'predicted_mrna', 'first_changed_aa', 'stop_codon_cdna_1based',
           'stop_to_last_junction_nt', 'nmd_prediction', 'protein_prediction_status', 'assumptions']


def read_table(output, name):
    with (Path(output) / (name + '.csv')).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def write_table(output, name, rows, columns):
    target = Path(output) / name
    temporary = target.with_suffix(target.suffix + '.partial')
    with temporary.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(target)


def model_loss(transcript, removed, fa, chrom):
    """Remove whole annotated exons, preserving strand and the annotated start.

    Never invent a downstream initiation site or a sequence at an unknown end.
    A terminal-loss peptide is only the translated known fragment, not a complete
    mutant protein. Internal joining remains conditional without junction reads.
    """
    seq, g2t, start, end, _ = rnako.template(transcript, fa, chrom)
    removed = set(removed)
    indices = {g2t[p] for a, b in removed for p in range(a, b) if p in g2t}
    remaining = [ex for ex in transcript.ordered_exons() if ex not in removed]
    new_seq = ''.join(base for i, base in enumerate(seq) if i not in indices)
    cds_positions = {g2t[p] for a, b in transcript.cds for p in range(a, b) if p in g2t}
    coding_deleted = len(cds_positions & indices)
    out = dict(effect='noncoding_region_loss', deleted_coding_nt=coding_deleted,
               coding_loss_fraction=round(coding_deleted / len(cds_positions), 5) if cds_positions else 0,
               reference_protein='', predicted_protein='', predicted_mrna=new_seq,
               first_changed_aa='', stop_codon_cdna_1based='', stop_to_last_junction_nt='',
               nmd_prediction='unknown', protein_prediction_status='not_available')
    if start is None or seq[start:start+3] != 'ATG':
        out['effect'] = 'coding_annotation_unresolved' if coding_deleted else out['effect']
        return out
    reference, _ = rnako.translate(seq[start:])
    out['reference_protein'] = reference
    if any(p in indices for p in range(start, start+3)):
        out.update(effect='start_region_lost', protein_prediction_status='start_lost; alternative_initiation_not_predicted')
        return out
    new_start = start - sum(p < start for p in indices)
    peptide, stop = rnako.translate(new_seq[new_start:])
    out['predicted_protein'] = peptide
    terminal = transcript.ordered_exons()[-1] in removed
    out['effect'] = ('noncoding_region_loss' if not coding_deleted else
                     'terminal_coding_region_loss' if terminal else
                     'region_loss_frameshift' if coding_deleted % 3 else 'region_loss_inframe')
    out['protein_prediction_status'] = ('conditional_known_fragment; unknown_3prime_end' if terminal else
                                       'conditional_exon_join_model; other_edits_not_phased')
    if reference != peptide:
        out['first_changed_aa'] = next((i+1 for i, (a, b) in enumerate(zip(reference, peptide)) if a != b),
                                       min(len(reference), len(peptide))+1)
    if stop is not None and not terminal:
        stop_pos = new_start + stop
        out['stop_codon_cdna_1based'] = stop_pos + 1
        boundaries = [sum(b-a for a, b in remaining[:i]) for i in range(1, len(remaining))]
        if boundaries:
            distance = boundaries[-1] - stop_pos - 3
            out['stop_to_last_junction_nt'] = distance
            out['nmd_prediction'] = ('possible_under_model' if distance > 55 else
                                     'unlikely_under_model' if distance <= 50 else 'borderline_under_model')
    return out


def assess_regions(samples, exons, min_depth=10):
    """Require WT exon coverage, valid libraries and retained KO coverage.

    A >=90% normalized depth reduction is usable even with a few residual reads.
    A >=80% reduction additionally requires <50% covered bases. Global RNA loss,
    absent libraries, and exons unexpressed in WT cannot generate a region call.
    """
    libraries = {s['srx_id']: float(s.get('unique_primary_library_fragments') or 0) for s in samples}
    def key(row):
        return (int(row['exon_start_1based'])-1, int(row['exon_end_1based']))
    wt = defaultdict(list)
    for row in exons:
        lib = libraries.get(row['srx_id'], 0)
        if (row['metadata_label'] == 'WT' and lib > 0
                and float(row['sufficient_depth_fraction']) >= .8
                and float(row['mean_depth']) >= min_depth):
            wt[key(row)].append(float(row['mean_depth']) / lib)
    loss, anchors = {}, {}
    for sample in samples:
        sid = sample['srx_id']
        if sample['metadata_label'] != 'KO' or libraries.get(sid, 0) <= 0:
            continue
        ratios, retained = {}, []
        for row in exons:
            interval = key(row)
            if row['srx_id'] != sid or not wt[interval]:
                continue
            ratio = float(row['mean_depth']) / libraries[sid] / mean(wt[interval])
            if (float(row['mean_depth']) >= min_depth and float(row['sufficient_depth_fraction']) >= .8
                    and ratio >= .5):
                retained.append(interval)
            if ratio <= .1 or (ratio <= .2 and float(row['covered_fraction']) < .5):
                ratios[interval] = (ratio, len(wt[interval]))
        if retained:
            loss[sid] = ratios
            anchors[sid] = retained
    return loss, anchors


def write_region_loss(output, fasta, gtf, gene, check_pause=None):
    output = Path(output)
    checkpoint = check_pause or (lambda: None)
    samples = read_table(output, 'samples')
    exons = read_table(output, 'exon_coverage')
    junctions = read_table(output, 'junctions')
    settings = json.loads((output / 'run.json').read_text(encoding='utf-8'))
    thresholds = settings.get('thresholds', {})
    losses, anchors = assess_regions(samples, exons, thresholds.get('min_event_depth', 10))
    transcripts = rnako.load_gene(gtf, gene)
    models = []
    with rnako.pysam.FastaFile(str(fasta)) as fa:
        # Prefer the longest CDS with a resolved canonical initiation site. Some
        # annotations have incomplete/non-ATG starts; keep their models marked
        # unresolved instead of allowing them to hide the resolved CDS models.
        templates = {t.tid: rnako.template(t, fa, rnako.resolve_reference_name(t.chrom, fa))
                     for t in transcripts.values()}
        resolved = {t.tid for t in transcripts.values() if templates[t.tid][2] is not None
                    and templates[t.tid][0][templates[t.tid][2]:templates[t.tid][2]+3] == 'ATG'}
        max_cds = max((sum(b-a for a, b in t.cds) for t in transcripts.values() if t.tid in resolved), default=0)
        for sid, lost in losses.items():
            for t in transcripts.values():
                checkpoint()
                ordered = t.ordered_exons()
                # Jointly model all depleted exons of this transcript; never create
                # inconsistent independent peptides for a multi-exon loss.
                removed = [ex for ex in ordered if ex in lost]
                if not removed or len(removed) == len(ordered):
                    continue
                chrom = rnako.resolve_reference_name(t.chrom, fa)
                prediction = model_loss(t, removed, fa, chrom)
                first, last = min(ordered.index(ex) for ex in removed), max(ordered.index(ex) for ex in removed)
                support = 0
                if first > 0 and last+1 < len(ordered) and removed == ordered[first:last+1]:
                    left, right = ordered[first-1], ordered[last+1]
                    junction = (min(left[1], right[1])+1, max(left[0], right[0]))
                    support = max((int(row['support_fragments']) for row in junctions
                        if row['srx_id'] == sid and (int(row['intron_start_1based']), int(row['intron_end_1based'])) == junction
                        and row['comparison_status'] == 'enriched_vs_comparators'), default=0)
                major = t.tid in resolved and max_cds > 0 and sum(b-a for a, b in t.cds) == max_cds
                disruptive = (prediction['effect'] == 'start_region_lost' or
                              prediction['deleted_coding_nt'] >= 30 and
                              (prediction['coding_loss_fraction'] >= .1 or prediction['effect'] == 'region_loss_frameshift'))
                ko_support = major and disruptive
                event = hashlib.sha256(f'{sid}:{t.tid}:{removed}'.encode()).hexdigest()[:16]
                models.append(dict(event_id=event, srx_id=sid, metadata_label='KO', transcript_id=t.tid,
                    chrom=t.chrom, strand=t.strand, lost_exons=';'.join(str(ordered.index(ex)+1) for ex in removed),
                    lost_regions=';'.join(f'{t.chrom}:{a+1}-{b}' for a, b in removed),
                    wt_reference_samples=min(lost[ex][1] for ex in removed),
                    maximum_normalized_depth_ratio=round(max(lost[ex][0] for ex in removed), 5),
                    retained_control_regions=';'.join(f'{t.chrom}:{a+1}-{b}' for a, b in anchors[sid]),
                    evidence='junction_supported' if support >= thresholds.get('min_junction', 3) else 'coverage_only',
                    junction_support=support, ko_support='yes' if ko_support else 'no', **prediction,
                    assumptions='RNA領域欠損の推定。DNA欠失・実際の接合境界・タンパク質の存在は未確定。'
                                '注釈エクソン境界で欠損を仮定した配列。末端欠損は既知領域の断片のみ。'))
    apply_region_classes(samples, exons, models, thresholds.get('min_event_depth', 10))
    for row in exons:
        interval = (int(row['exon_start_1based'])-1, int(row['exon_end_1based']))
        match = losses.get(row['srx_id'], {}).get(interval)
        row['region_loss_status'] = 'suspected_region_loss' if match else ''
        row['region_normalized_depth_ratio'] = round(match[0], 5) if match else ''
        row['region_wt_reference_samples'] = match[1] if match else ''
    write_table(output, 'exon_coverage.csv', exons, list(exons[0]) if exons else
                ['srx_id','metadata_label','exon_start_1based','exon_end_1based'])
    write_table(output, 'region_deletions.csv', models, COLUMNS)
    write_table(output, 'samples.csv', samples, list(samples[0]))
    for filename, column in [('region_proteins.fasta', 'predicted_protein'), ('region_transcripts.fasta', 'predicted_mrna')]:
        target = output / filename
        temporary = target.with_suffix('.fasta.partial')
        with temporary.open('w', encoding='utf-8') as handle:
            for row in models:
                if row[column]:
                    handle.write(f'>{row["event_id"]}|{row["srx_id"]}|{row["transcript_id"]}|{row["evidence"]}|{row["protein_prediction_status"]}\n')
                    sequence = row[column]
                    handle.write('\n'.join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + '\n')
        temporary.replace(target)
    settings['region_loss_model'] = {'version': 1, 'normalized_depth_ratio': .1,
        'partial_coverage_ratio': .2, 'retained_normalized_ratio': .5,
        'minimum_coding_loss_nt': 30, 'minimum_coding_loss_fraction': .1,
        'classification': 'RNA KO_like; longest canonical-start annotated CDS affected; DNA deletion unconfirmed'}
    target = output / 'run.json'
    temporary = target.with_suffix('.json.partial')
    temporary.write_text(json.dumps(settings, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    temporary.replace(target)
    from structure import write_structure
    write_structure(output, fasta, gtf, gene, transcripts)
    return models
