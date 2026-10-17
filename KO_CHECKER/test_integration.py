"""Indexed BAM regression tests for sample-specific splice classification."""
import csv
from pathlib import Path
import tempfile
import unittest

import rnako
from pipeline import write_exon_coverage, atomic_json, PipelinePaused


@unittest.skipIf(rnako.pysam is None, 'pysam required')
class BamIntegrationTests(unittest.TestCase):
    def test_splice_classification_is_independent_of_sample_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sequence = 'ATGAAACCC' + 'N' * 11 + 'CCCGGG' + 'N' * 14 + 'TTTTAACCC'
            fasta = root / 'reference.fa'
            fasta.write_text('>chr1\n' + sequence + '\n')
            gtf = root / 'reference.gtf'
            attributes = 'gene_id "G1"; gene_name "Test"; transcript_id "T1"; gene_type "protein_coding";'
            with gtf.open('w') as handle:
                for start, end in ((1, 9), (21, 26), (41, 49)):
                    handle.write(f'chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attributes}\n')
                    handle.write(f'chr1\ttest\tCDS\t{start}\t{min(end, 43)}\t.\t+\t0\t{attributes}\n')
            for sid, cigar, bases in (
                ('WT', '9M11N6M14N9M', sequence[:9] + sequence[20:26] + sequence[40:]),
                ('KO', '9M31N9M', sequence[:9] + sequence[40:]),
                ('WT_rare', '9M11N6M14N9M', sequence[:9] + sequence[20:26] + sequence[40:]),
            ):
                path = root / f'{sid}.bam'
                with rnako.pysam.AlignmentFile(str(path), 'wb', header={
                    'HD': {'VN': '1.6', 'SO': 'coordinate'},
                    'SQ': [{'SN': '1', 'LN': len(sequence)}],
                }) as bam:
                    for i in range(12 if sid == "KO" else 104 if sid == "WT_rare" else 100):
                        read = rnako.pysam.AlignedSegment()
                        read.query_name = f'{sid}_{i}'
                        read.query_sequence = sequence[:9] + sequence[40:] if sid == "WT_rare" and i < 4 else bases
                        read.reference_id = 0
                        read.reference_start = 0
                        read.mapping_quality = 60
                        read.cigarstring = "9M31N9M" if sid == "WT_rare" and i < 4 else cigar
                        read.query_qualities = rnako.pysam.qualitystring_to_array('I' * len(read.query_sequence))
                        read.set_tag('NH', 1)
                        bam.write(read)
                rnako.pysam.index(str(path))
            for order in (('WT', 'WT_rare', 'KO'), ('KO', 'WT_rare', 'WT')):
                with self.subTest(order=order):
                    samples = root / 'input.csv'
                    with samples.open('w', newline='') as handle:
                        writer = csv.writer(handle)
                        writer.writerow(['srx_id', 'metadata_label', 'bam_path'])
                        for sid in order:
                            writer.writerow([sid, 'WT' if sid.startswith('WT') else 'KO', root / f'{sid}.bam'])
                    out = root / ('out_' + '_'.join(order))
                    args = rnako.parser().parse_args([
                        '--gene', 'Test', '--samples', str(samples),
                        '--fasta', str(fasta), '--gtf', str(gtf),
                        '--out', str(out), '--min-gene-fragments', '1',
                    ])
                    rnako.run(args)
                    with (out / 'samples.csv').open() as handle:
                        rows = {row['srx_id']: row for row in csv.DictReader(handle)}
                    self.assertEqual(rows['KO']['predicted_class'], 'KO_like')
                    self.assertEqual(rows['WT']['predicted_class'], 'WT_like')
                    self.assertEqual(rows['WT_rare']['predicted_class'], 'WT_like')
                    self.assertEqual(float(rows['WT']['expression_relative_to_wt']), 1.0)
                    with samples.open() as handle:
                        coverage_samples = list(csv.DictReader(handle))
                    control = next(r for r in coverage_samples if r['srx_id'] == 'WT')
                    coverage_samples.append({**control, 'srx_id': 'CONTROL', 'metadata_label': 'KO'})
                    write_exon_coverage(coverage_samples, out, fasta, gtf, 'Test')
                    with (out / 'exon_coverage.csv').open() as handle:
                        coverage = list(csv.DictReader(handle))
                    middle = next(r for r in coverage if r['srx_id'] == 'KO' and r['exon_start_1based'] == '21')
                    self.assertEqual(middle['comparison_status'], 'suspected_local_depletion')
                    self.assertEqual(middle['coverage_status'], 'no_coverage')
                    self.assertEqual(float(middle['relative_depth_vs_wt']), 0)
                    self.assertEqual(middle['skipping_junction_fragments'], '12')
                    self.assertEqual(middle['longest_zero_depth_region'], 'chr1:21-26')
                    outer = next(r for r in coverage if r['srx_id'] == 'KO' and r['exon_start_1based'] == '1')
                    self.assertEqual(outer['comparison_status'], 'insufficient_sample_coverage')
                    self.assertEqual(outer['mean_depth'], '12.0')
                    self.assertTrue(all(r['comparison_status'] == 'no_specific_depletion'
                                        for r in coverage if r['srx_id'] == 'CONTROL'))
                    previous_report = (out / 'exon_coverage.csv').read_bytes()
                    checkpoints = 0
                    def pause_during_coverage():
                        nonlocal checkpoints
                        checkpoints += 1
                        if checkpoints == 3:
                            raise PipelinePaused('pause during pileup')
                    with self.assertRaises(PipelinePaused):
                        write_exon_coverage(coverage_samples, out, fasta, gtf, 'Test', pause_during_coverage)
                    self.assertEqual((out / 'exon_coverage.csv').read_bytes(), previous_report)
                    settings = {'thresholds': {'min_mapq': 0, 'min_bq': 20, 'min_event_depth': 1000}}
                    atomic_json(out / 'run.json', settings)
                    write_exon_coverage(coverage_samples, out, fasta, gtf, 'Test')
                    with (out / 'exon_coverage.csv').open() as handle:
                        insufficient = [r for r in csv.DictReader(handle) if r['metadata_label'] == 'KO']
                    self.assertTrue(all(r['comparison_status'] == 'insufficient_wt_coverage' for r in insufficient))
                    self.assertTrue(all(r['relative_depth_vs_wt'] == '' for r in insufficient))


if __name__ == '__main__':
    unittest.main()
