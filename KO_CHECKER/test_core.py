import csv
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

from rnako import Transcript, junction_effect, load_gene, protein_effect, rc, read_samples, template, resolve_contig


class FakeReference:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, chrom, start, end):
        return self.sequence[start:end]


class ConsequenceTests(unittest.TestCase):
    def setUp(self):
        # Two exons, first 9 nt ATG AAA CCC, second 9 nt GGG TAA AAA.
        bases = list("N" * 60)
        bases[0:9] = list("ATGAAACCC")
        bases[20:29] = list("GGGTAAAAA")
        self.fa = FakeReference("".join(bases))
        self.tx = Transcript("T1", "G1", "Test", "chr1", "+", "protein_coding", [(0, 9), (20, 29)], [(0, 9), (20, 26)])

    def test_synonymous_variant(self):
        c = protein_effect(self.tx, {"pos": 6, "ref": "A", "alt": "G"}, self.fa, "chr1")
        self.assertEqual(c["effect"], "synonymous_or_same_predicted_protein")
        self.assertEqual(c["reference_protein"], "MKPG")
        self.assertEqual(c["predicted_protein"], "MKPG")

    def test_frameshift_insertion(self):
        c = protein_effect(self.tx, {"pos": 5, "ref": "A", "alt": "AT"}, self.fa, "chr1")
        self.assertEqual(c["effect"], "frameshift")
        self.assertNotEqual(c["reference_protein"], c["predicted_protein"])

    def test_stop_gained(self):
        # AAA -> TAA: a genuine premature stop, not an in-frame deletion.
        c = protein_effect(self.tx, {"pos": 4, "ref": "A", "alt": "T"}, self.fa, "chr1")
        self.assertEqual(c["effect"], "stop_gained")
        self.assertEqual(c["predicted_protein"], "M")

    def test_inframe_deletion_is_not_stop_gain(self):
        c = protein_effect(self.tx, {"pos": 3, "ref": "GAAA", "alt": "G"}, self.fa, "chr1")
        self.assertEqual(c["effect"], "inframe_deletion")

    def test_exon_skip_changes_protein(self):
        # Skipping a terminal exon containing stop may lead to no predicted stop.
        bases = list("N"*80)
        bases[0:6] = list("ATGAAA")
        bases[20:23] = list("CCC")
        bases[40:46] = list("GGGTAA")
        fa = FakeReference("".join(bases))
        tx = Transcript("T2", "G1", "Test", "chr1", "+", "protein_coding", [(0,6),(20,23),(40,46)], [(0,6),(20,23),(40,43)])
        c = junction_effect(tx, (6,40), fa, "chr1")
        self.assertEqual(c["effect"], "exon_skip_inframe")
        self.assertEqual(c["predicted_protein"], "MKG")

    def test_reverse_strand_substitution(self):
        # Genomic CAT reverse complements to ATG (start), then TTT -> AAA, TTA -> TAA.
        bases = list("N"*40)
        bases[20:29] = list("TTATTT CAT".replace(" ", ""))
        fa = FakeReference("".join(bases))
        tx = Transcript("Tminus", "G1", "Test", "chr1", "-", "protein_coding", [(20,29)], [(20,26), (26,29)])
        seq, _, start, _, _ = template(tx, fa, "chr1")
        self.assertEqual(seq[start:start+3], "ATG")

    def test_rc(self):
        self.assertEqual(rc("ACGT"), "ACGT")


class InputTests(unittest.TestCase):
    def test_gene_symbol_case_fallback_preserves_exact_match_and_ambiguity(self):
        with tempfile.TemporaryDirectory() as temp:
            gtf = Path(temp) / 'test.gtf'
            line = 'chr1\tx\texon\t1\t9\t.\t+\t.\tgene_id "G1"; gene_name "Zfp445"; transcript_id "T1";\n'
            gtf.write_text(line)
            for query in ('Zfp445', 'ZFP445', 'zfp445', 'G1'):
                with self.subTest(query=query):
                    self.assertEqual(load_gene(gtf, query)['T1'].gene_name, 'Zfp445')
            with self.assertRaisesRegex(ValueError, 'not found'):
                load_gene(gtf, 'ZFP444')
            gtf.write_text(line + line.replace('G1', 'G2').replace('T1', 'T2').replace('Zfp445', 'ZFP445'))
            self.assertEqual(set(load_gene(gtf, 'Zfp445')), {'T1'})
            self.assertEqual(set(load_gene(gtf, 'ZFP445')), {'T2'})
            with self.assertRaisesRegex(ValueError, 'ambiguous'):
                load_gene(gtf, 'zfp445')

    def test_contig_alias_keeps_length_check(self):
        bam = SimpleNamespace(references=("1",), get_reference_length=lambda name: 100)
        fasta = SimpleNamespace(references=("chr1",), get_reference_length=lambda name: 100)
        self.assertEqual(resolve_contig("chr1", bam, fasta), "1")
        fasta.get_reference_length = lambda name: 101
        with self.assertRaisesRegex(ValueError, "length differs"):
            resolve_contig("chr1", bam, fasta)


    def test_gtf_and_explicit_path(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            gtf = temp / "test.gtf"
            gtf.write_text('chr1\tx\texon\t1\t9\t.\t+\t.\tgene_id "G1"; gene_name "Test"; transcript_id "T1"; gene_type "protein_coding";\n'
                           'chr1\tx\tCDS\t1\t9\t.\t+\t0\tgene_id "G1"; gene_name "Test"; transcript_id "T1";\n')
            tx = load_gene(gtf, "Test")
            self.assertEqual(tx["T1"].exons, [(0, 9)])
            (temp / "SRX1.bam").touch()
            samples = temp / "samples.csv"
            with samples.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["srx_id", "metadata_label", "bam_path"])
                w.writerow(["SRX1", "KO", str(temp / "SRX1.bam")])
            args = SimpleNamespace(samples=str(samples), bam_pattern="unused", wt="", ko="", kd="", oe="")
            self.assertEqual(read_samples(args)[0]["metadata_label"], "KO")


if __name__ == "__main__":
    unittest.main()
