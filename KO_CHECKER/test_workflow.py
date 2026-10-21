"""Download/preprocessing cache regressions and offline real-tool mapping tests."""
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import rnako
from pipeline import Pipeline, PipelinePaused, atomic_json, digest, mark_complete


class PreprocessingTests(unittest.TestCase):
    def engine(self, root, threads=2):
        return Pipeline(root / 'data', root / 'references', threads, lambda **kw: None, root / 'pipeline.log')

    @staticmethod
    def fake_fastp(args, stdout=None):
        for flag in ('-o', '-O'):
            if flag in args:
                with gzip.open(args[args.index(flag) + 1], 'wt') as handle:
                    handle.write('@r\n' + 'ACGT' * 10 + '\n+\n' + 'I' * 40 + '\n')
        Path(args[args.index('-h') + 1]).write_text('<html>QC</html>')
        atomic_json(args[args.index('-j') + 1], {'summary': {'after_filtering': {'total_reads': 2}}})

    def raw(self, root, paired=True):
        paths = [root / ('read' + suffix + '.fastq') for suffix in (('_1', '_2') if paired else ('',))]
        for path in paths:
            path.write_text('@r\n' + 'ACGT' * 10 + '\n+\n' + 'I' * 40 + '\n')
        return paths

    def test_old_raw_fastq_is_reused_without_downloading_again(self):
        with tempfile.TemporaryDirectory() as folder:
            engine = self.engine(Path(folder))
            directory = engine.data / 'fastq' / 'SRR1'
            directory.mkdir(parents=True)
            raw = directory / 'SRR1.fastq'
            raw.write_text('@r\nACGT\n+\nIIII\n')
            mark_complete(directory, {'version': 1, 'run': 'SRR1', 'layout': 'SINGLE',
                                      'fasterq': 'fq1', 'mode': 'split-files'}, [raw])
            with patch.object(engine, 'command') as command:
                result = engine.fastqs({'run_accession': 'SRR1', 'library_layout': 'SINGLE'}, {'fasterq-dump': 'fq1'})
            command.assert_not_called()
            self.assertEqual(result, [raw])

    def test_fastp_handles_both_layouts_and_reuses_complete_results(self):
        for paired in (True, False):
            with self.subTest(paired=paired), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                engine = self.engine(root, 23)
                raw = self.raw(root, paired)
                row = {'run_accession': 'SRR1', 'library_layout': 'PAIRED' if paired else 'SINGLE'}
                with patch.object(engine, 'command', side_effect=self.fake_fastp) as command:
                    result = engine.trimmed_fastqs(row, raw, {'fastp': 'v1'})
                    args = command.call_args.args[0]
                    self.assertEqual('--detect_adapter_for_pe' in args, paired)
                    self.assertEqual('-I' in args, paired)
                    self.assertEqual(args[args.index('-w') + 1], 16)
                    engine.threads = 4
                    self.assertEqual(engine.trimmed_fastqs(row, raw, {'fastp': 'v1'}), result)
                    self.assertEqual(command.call_count, 1)
                    self.assertEqual(engine.preprocessing['SRR1']['fastp_threads'], 16)
                    raw[0].write_text(raw[0].read_text() + '@new\nACGT\n+\nIIII\n')
                    changed = engine.trimmed_fastqs(row, raw, {'fastp': 'v1'})
                    self.assertNotEqual(changed, result)
                    self.assertEqual(command.call_count, 2)
                    upgraded = engine.trimmed_fastqs(row, raw, {'fastp': 'v2'})
                    self.assertNotEqual(upgraded, changed)

    def test_missing_qc_or_mate_is_never_marked_complete(self):
        for missing in ('-O', '-h', '-j'):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                engine = self.engine(root)
                raw = self.raw(root)
                def incomplete(args, stdout=None):
                    self.fake_fastp(args)
                    Path(args[args.index(missing) + 1]).unlink()
                with patch.object(engine, 'command', side_effect=incomplete), self.assertRaisesRegex(RuntimeError, 'fastp'):
                    engine.trimmed_fastqs({'run_accession': 'SRR1', 'library_layout': 'PAIRED'}, raw, {'fastp': 'v1'})
                self.assertEqual(list(engine.data.rglob('complete.json')), [])

    def test_all_reads_filtered_produces_actionable_error(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            engine = self.engine(root)
            def filtered(args, stdout=None):
                self.fake_fastp(args)
                atomic_json(args[args.index('-j') + 1], {'summary': {'after_filtering': {'total_reads': 0}}})
            with patch.object(engine, 'command', side_effect=filtered), self.assertRaisesRegex(RuntimeError, '有効なリード'):
                engine.trimmed_fastqs({'run_accession': 'SRR1', 'library_layout': 'SINGLE'}, self.raw(root, False), {'fastp': 'v1'})
            self.assertEqual(list(engine.data.rglob('complete.json')), [])

    def test_pause_does_not_publish_trimmed_outputs_and_resume_retries(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            engine = self.engine(root)
            raw = self.raw(root, False)
            row = {'run_accession': 'SRR1', 'library_layout': 'SINGLE'}
            def pause(args, stdout=None):
                self.fake_fastp(args)
                raise PipelinePaused('pause fastp')
            with patch.object(engine, 'command', side_effect=pause), self.assertRaises(PipelinePaused):
                engine.trimmed_fastqs(row, raw, {'fastp': 'v1'})
            self.assertEqual(list(engine.data.rglob('complete.json')), [])
            self.assertEqual(list(engine.data.rglob('*.trim.fastq.gz')), [])
            with patch.object(engine, 'command', side_effect=self.fake_fastp):
                paths = engine.trimmed_fastqs(row, raw, {'fastp': 'v1'})
            self.assertTrue(paths[0].is_file())

    def test_old_untrimmed_bam_is_retained_but_never_reused(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            engine = self.engine(root)
            index = root / 'index' / 'genome'
            atomic_json(index.parent / 'complete.json', {'reference': 'mm10'})
            old = {'version': 1, 'run': 'SRR1', 'layout': 'SINGLE', 'index': {'reference': 'mm10'},
                   'hisat2': 'h1', 'samtools': 's1', 'options': ['--dta', '--rg-id', 'SRR1']}
            directory = engine.data / 'bam' / 'SRR1' / digest(old)
            directory.mkdir(parents=True)
            bam = directory / 'SRR1.bam'
            bam.write_bytes(b'original')
            mark_complete(directory, old, [bam])
            with patch.object(engine, 'fastqs', side_effect=RuntimeError('new preprocessing required')) as fastqs, \
                 self.assertRaisesRegex(RuntimeError, 'new preprocessing'):
                engine.mapped_run({'run_accession': 'SRR1', 'library_layout': 'SINGLE'}, index,
                    {'hisat2': 'h1', 'samtools': 's1', 'fastp': 'p1', 'fasterq-dump': 'q1'})
            fastqs.assert_called_once()
            self.assertEqual(bam.read_bytes(), b'original')


@unittest.skipUnless(rnako.pysam is not None and all(shutil.which(n) for n in ('fastp', 'hisat2', 'hisat2-build', 'samtools')),
                     'fastp/HISAT2/samtools and pysam required')
class RealToolTests(unittest.TestCase):
    def test_real_fastp_mapping_index_and_bam_cache_for_both_layouts(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            rng = random.Random(51)
            genome = ''.join(rng.choice('ACGT') for _ in range(2000))
            genome = genome[:900] + genome[100:200] + genome[1000:]
            fasta = root / 'reference.fa'
            fasta.write_text('>chr1\n' + genome + '\n')
            rnako.pysam.faidx(str(fasta))
            engine = Pipeline(root / 'data', root / 'references', 2, lambda **kw: None, root / 'pipeline.log')
            versions = {name: subprocess.check_output([name, '--version'], stderr=subprocess.STDOUT).decode(errors='replace')
                        for name in ('fastp', 'fasterq-dump', 'hisat2', 'hisat2-build', 'samtools')}
            index = engine.index(fasta, versions)
            for layout in ('SINGLE', 'PAIRED'):
                row = {'run_accession': 'SRR1' if layout == 'SINGLE' else 'SRR2', 'library_layout': layout}
                raw = [root / f'{layout}_1.fastq']
                if layout == 'PAIRED':
                    raw.append(root / f'{layout}_2.fastq')
                for mate, path in enumerate(raw):
                    bases = genome[300:375] if mate == 0 else rnako.rc(genome[450:525])
                    with path.open('w') as handle:
                        for i in range(20):
                            handle.write(f'@read{i}\n{bases}\n+\n' + 'I' * len(bases) + '\n')
                with self.subTest(layout=layout), patch.object(engine, 'fastqs', return_value=raw):
                    bam = engine.mapped_run(row, index, versions)
                with rnako.pysam.AlignmentFile(str(bam), 'rb') as handle:
                    self.assertTrue(handle.has_index())
                    reads = list(handle.fetch())
                    self.assertTrue(reads)
                    self.assertTrue(all(r.get_tag('NH') == 1 for r in reads))
                    if layout == 'PAIRED':
                        self.assertTrue(all(r.is_proper_pair for r in reads))
                self.assertTrue(Path(engine.preprocessing[row['run_accession']]['fastp_json']).is_file())
                with patch.object(engine, 'fastqs') as raw_again, patch.object(engine, 'trimmed_fastqs') as trim_again:
                    self.assertEqual(engine.mapped_run(row, index, versions), bam)
                    raw_again.assert_not_called()
                    trim_again.assert_not_called()
            repeated = root / 'repeated.fastq'
            with repeated.open('w') as handle:
                for i in range(20):
                    handle.write(f'@repeat{i}\n{genome[100:175]}\n+\n' + 'I' * 75 + '\n')
            with patch.object(engine, 'fastqs', return_value=[repeated]):
                repeated_bam = engine.mapped_run({'run_accession': 'SRR3', 'library_layout': 'SINGLE'}, index, versions)
            with rnako.pysam.AlignmentFile(str(repeated_bam), 'rb') as handle:
                repeat_reads = list(handle.fetch())
                self.assertTrue(any(r.is_secondary for r in repeat_reads))
                self.assertTrue(all(r.get_tag('NH') >= 2 for r in repeat_reads))
                self.assertFalse(any(rnako.read_ok(r, 0) for r in repeat_reads))
            log = engine.log.read_text()
            self.assertIn('-k 100 --secondary', log)
            self.assertIn('--no-mixed --no-discordant', log)
            self.assertNotIn('--dta', log)
            self.assertIn('.trim.fastq.gz', log)
            wrong = root / 'wrong.bam'
            with rnako.pysam.AlignmentFile(str(wrong), 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 1999}]}):
                pass
            with self.assertRaisesRegex(RuntimeError, '全contig長'):
                engine.check_bam(wrong)


if __name__ == '__main__':
    unittest.main()
