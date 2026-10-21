"""Regression tests for scientific claims in conclusion-first reports."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

from pipeline import validate_request
from reporting import build_conclusion, permutation_test, read_conclusion, write_conclusion


def sample(sid, label, cpm, verdict='undetermined'):
    return {'srx_id': sid, 'metadata_label': label, 'gene': 'Test',
            'gene_unique_fragment_cpm': cpm, 'unique_primary_library_fragments': 1000,
            'predicted_class': verdict, 'coverage_warning': ''}


class ConclusionTests(unittest.TestCase):
    def setUp(self):
        self.tables = {'samples': [sample('W', 'WT', 100), sample('K', 'KO', 20)]}

    def test_expression_is_compared_when_all_sequence_calls_are_undetermined(self):
        report = build_conclusion(self.tables)
        self.assertEqual(report['sequence']['status'], 'undetermined')
        self.assertEqual(report['expression']['ratio'], .2)
        self.assertEqual(report['expression']['decrease_percent'], 80)
        self.assertIn('80.0%低下', report['expression']['text'])
        self.assertEqual(report['expression']['test']['status'], 'insufficient_replicates')
        self.assertEqual(report['nmd']['status'], 'unresolved')

    def test_sequence_coverage_does_not_prevent_expression_comparison(self):
        self.tables['samples'][1]['coverage_warning'] = 'low_exonic_fragments'
        report = build_conclusion(self.tables)
        self.assertEqual(report['sequence']['status'], 'insufficient_coverage')
        self.assertEqual(report['expression']['decrease_percent'], 80)

    def test_all_labelled_wt_samples_form_the_baseline(self):
        self.tables['samples'].append(sample('W2', 'WT', 300, 'WT_like'))
        report = build_conclusion(self.tables)
        self.assertEqual(report['expression']['wt_mean_cpm'], 200)
        self.assertEqual(report['expression']['ratio'], .1)

    def test_zero_ko_expression_is_a_valid_complete_reduction(self):
        self.tables['samples'][1]['gene_unique_fragment_cpm'] = 0
        report = build_conclusion(self.tables)
        self.assertEqual(report['expression']['decrease_percent'], 100)
        self.assertEqual(report['expression']['ratio'], 0)

    def test_zero_wt_and_missing_library_cannot_produce_a_percentage(self):
        self.tables['samples'][0]['gene_unique_fragment_cpm'] = 0
        self.assertIsNone(build_conclusion(self.tables)['expression']['ratio'])
        self.tables['samples'][0]['unique_primary_library_fragments'] = 0
        report = build_conclusion(self.tables)
        self.assertIsNone(report['expression']['ratio'])
        self.assertEqual(report['expression']['excluded_samples'], ['W'])

    def test_increases_are_not_described_as_reductions(self):
        self.tables['samples'][1]['gene_unique_fragment_cpm'] = 200
        report = build_conclusion(self.tables)
        self.assertIn('100.0%増加', report['expression']['text'])
        self.assertNotEqual(report['nmd']['status'], 'consistent')

    def test_single_replicates_do_not_get_p_values_even_with_confirmation(self):
        report = build_conclusion(self.tables, True)
        self.assertEqual(report['expression']['test']['status'], 'insufficient_replicates')
        self.assertNotIn('p_value', report['expression']['test'])

    def test_large_read_counts_are_not_treated_as_biological_replicates(self):
        self.tables['samples'][0]['unique_primary_library_fragments'] = 100000000
        report = build_conclusion(self.tables, True)
        self.assertEqual(report['expression']['test']['status'], 'insufficient_replicates')

    def test_exact_permutation_reference_case_and_identical_groups(self):
        result = permutation_test([100] * 4, [10] * 4)
        self.assertAlmostEqual(result['p_value'], 2 / 70)
        self.assertTrue(result['significant'])
        self.assertEqual(permutation_test([10] * 4, [10] * 4)['p_value'], 1)
        self.assertFalse(permutation_test([100] * 2, [0] * 2)['significant'])

    def test_randomized_test_is_nonzero_and_independent_of_sample_order(self):
        wt, ko = list(range(100, 120)), list(range(20))
        result = permutation_test(wt, ko)
        self.assertEqual(result['method'], 'monte_carlo')
        self.assertGreater(result['p_value'], 0)
        self.assertEqual(result, permutation_test(list(reversed(wt)), list(reversed(ko))))

    def test_confirmed_groups_are_tested_but_unconfirmed_groups_are_not(self):
        self.tables['samples'] = ([sample(f'W{i}', 'WT', 100) for i in range(4)]
                                  + [sample(f'K{i}', 'KO', 10) for i in range(4)])
        self.assertEqual(build_conclusion(self.tables)['expression']['test']['status'], 'not_confirmed')
        test = build_conclusion(self.tables, True)['expression']['test']
        self.assertTrue(test['significant'])
        self.assertIn('探索的', test['text'])

    def test_nmd_requires_enriched_disruptive_sequence_and_reduction_in_same_sample(self):
        self.tables['events'] = [{'srx_id': 'K', 'event_id': 'E', 'comparison_status': 'enriched_vs_comparators'}]
        self.tables['consequences'] = [{'srx_id': 'K', 'event_id': 'E', 'effect': 'frameshift', 'nmd_prediction': 'possible'}]
        self.assertEqual(build_conclusion(self.tables)['nmd']['status'], 'consistent')
        self.tables['events'][0]['comparison_status'] = 'shared_or_not_enriched'
        self.assertEqual(build_conclusion(self.tables)['nmd']['status'], 'unresolved')
        self.tables['events'][0]['comparison_status'] = 'enriched_vs_comparators'
        self.tables['consequences'][0]['effect'] = 'synonymous_or_same_predicted_protein'
        self.assertEqual(build_conclusion(self.tables)['nmd']['status'], 'unresolved')

    def test_nmd_does_not_combine_sequence_in_one_sample_with_loss_in_another(self):
        self.tables['samples'].append(sample('K2', 'KO', 100))
        self.tables['events'] = [{'srx_id': 'K2', 'event_id': 'E', 'comparison_status': 'enriched_vs_comparators'}]
        self.tables['consequences'] = [{'srx_id': 'K2', 'event_id': 'E', 'effect': 'frameshift', 'nmd_prediction': 'possible'}]
        self.assertEqual(build_conclusion(self.tables)['nmd']['status'], 'unresolved')

    def test_exon_highlights_and_missing_exon_coverage(self):
        self.tables['exon_coverage'] = [{'srx_id': 'K', 'metadata_label': 'KO',
            'chrom': 'chr1', 'exon_start_1based': '21', 'exon_end_1based': '30',
            'relative_depth_vs_wt': '0.1', 'comparison_status': 'suspected_local_depletion'}]
        report = build_conclusion(self.tables)
        self.assertEqual(report['exons']['depleted_count'], 1)
        self.assertEqual(report['exons']['strongest'][0]['exon_start_1based'], '21')
        self.assertIn('統計検定による判定ではありません', report['exons']['text'])
        self.tables.pop('exon_coverage')
        self.assertIn('データがありません', build_conclusion(self.tables)['exons']['text'])

    def test_old_results_get_reports_without_rewriting_the_csv(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'samples.csv'
            with path.open('w', newline='', encoding='utf-8') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(self.tables['samples'][0]))
                writer.writeheader()
                writer.writerows(self.tables['samples'])
            original = path.read_bytes()
            report = read_conclusion(folder)
            self.assertEqual(report['expression']['decrease_percent'], 80)
            write_conclusion(folder)
            self.assertEqual(path.read_bytes(), original)
            self.assertIn('80.0%低下', (Path(folder) / 'conclusion.txt').read_text(encoding='utf-8'))
            self.assertEqual(json.loads((Path(folder) / 'conclusion.json').read_text(encoding='utf-8'))['expression']['ratio'], .2)

    def test_replicate_confirmation_must_be_an_explicit_boolean(self):
        for value in (True, False, 'true', 1):
            request = validate_request({'gene': 'Test', 'wt': 'SRX1', 'ko': 'SRX2', 'confirm_replicates': value})
            self.assertEqual(request.get('confirm_replicates', False), value is True)

    def test_conclusions_include_exons_after_the_browser_row_limit(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'exon_coverage.csv'
            with path.open('w', newline='', encoding='utf-8') as handle:
                writer = csv.writer(handle)
                writer.writerow(['srx_id', 'metadata_label', 'comparison_status', 'relative_depth_vs_wt',
                                 'chrom', 'exon_start_1based', 'exon_end_1based'])
                for _ in range(1000):
                    writer.writerow(['K', 'KO', 'no_specific_depletion', 1, 'chr1', 1, 10])
                writer.writerow(['K', 'KO', 'suspected_local_depletion', .1, 'chr1', 21, 30])
            report = write_conclusion(folder)
            self.assertEqual(report['exons']['depleted_count'], 1)
            self.assertIn('chr1:21-30', (Path(folder) / 'conclusion.txt').read_text(encoding='utf-8'))


if __name__ == '__main__':
    unittest.main()
