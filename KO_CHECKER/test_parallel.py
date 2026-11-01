"""RAM work area, chunked mapping and analyses running at the same time."""
import gzip
import json
import os
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

import pipeline
from pipeline import GIB, Pipeline, PipelinePaused, ReadChunks, Shared, atomic_json
from webapp import Jobs
from test_batch import CSV, completed


def fastq(path, names):
    with gzip.open(path, 'wt') as handle:
        for name in names:
            handle.write(f'@{name}\nACGT\n+\nIIII\n')


def chunk_names(path):
    return [line[1:].strip() for i, line in enumerate(path.read_text().splitlines(True)) if i % 4 == 0]


class ReadChunkTests(unittest.TestCase):
    def test_paired_blocks_keep_mates_together(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            names = [f'r{i}' for i in range(25)]
            fastq(root / 'a_1.fq.gz', names)
            fastq(root / 'a_2.fq.gz', names)
            chunks = ReadChunks([root / 'a_1.fq.gz', root / 'a_2.fq.gz'], root, size=10)
            chunks.BLOCK_LINES = 8
            blocks = []
            try:
                while (paths := chunks.next()) is not None:
                    first, second = (chunk_names(p) for p in paths)
                    self.assertEqual(first, second)
                    blocks.append(first)
            finally:
                chunks.close()
            self.assertEqual([len(b) for b in blocks], [10, 10, 5])
            self.assertEqual(sum(blocks, []), names)
            self.assertFalse(any(p.exists() for p in chunks.sets[0] + chunks.sets[1]))

    def test_mate_count_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            fastq(root / 'a_1.fq.gz', ['r1', 'r2', 'r3'])
            fastq(root / 'a_2.fq.gz', ['r1', 'r2'])
            chunks = ReadChunks([root / 'a_1.fq.gz', root / 'a_2.fq.gz'], root, size=10)
            try:
                with self.assertRaisesRegex(RuntimeError, '行数'):
                    chunks.next()
            finally:
                chunks.close()


class SharedResourceTests(unittest.TestCase):
    def engine(self, temp, shared, **kw):
        return Pipeline(temp, temp, 1, lambda **k: None, Path(temp) / 'log.txt', shared=shared, **kw)

    def test_second_analysis_waits_for_the_same_run_and_can_pause_while_waiting(self):
        with tempfile.TemporaryDirectory() as temp:
            shared = Shared(2)
            first = self.engine(temp, shared)
            entered, release, order = threading.Event(), threading.Event(), []

            def hold():
                with first.exclusive('run:SRR1', 'wait'):
                    order.append('first')
                    entered.set()
                    release.wait(5)
            thread = threading.Thread(target=hold)
            thread.start()
            self.assertTrue(entered.wait(5))
            stop = threading.Event()
            paused = self.engine(temp, shared, should_pause=stop.is_set)
            stop.set()
            with self.assertRaises(PipelinePaused):
                with paused.exclusive('run:SRR1', '待機'):
                    pass
            second = self.engine(temp, shared)
            waiter = threading.Thread(target=lambda: second.exclusive('run:SRR1', '待機').__enter__() or order.append('second'))
            waiter.start()
            time.sleep(0.2)
            self.assertEqual(order, ['first'])
            release.set()
            thread.join(5)
            waiter.join(5)
            self.assertEqual(order, ['first', 'second'])
            self.assertIn('待機', (Path(temp) / 'log.txt').read_text())

    def test_ram_reservations_are_shared_between_analyses(self):
        with tempfile.TemporaryDirectory() as temp:
            shared = Shared(2)
            usage = os.statvfs_result if False else None
            with patch('pipeline.shutil.disk_usage', return_value=type('U', (), {'total': 100 * GIB, 'free': 100 * GIB})()):
                a, b = self.engine(temp, shared, ram=temp), self.engine(temp, shared, ram=temp)
                a.reserve_ram('SRR1', 60 * GIB)
                with self.assertRaisesRegex(RuntimeError, '収まりません'):
                    b.reserve_ram('SRR2', 95 * GIB)
                done = threading.Event()
                threading.Thread(target=lambda: (b.reserve_ram('SRR2', 40 * GIB), done.set())).start()
                self.assertFalse(done.wait(0.5))
                a.shrink_ram(20 * GIB)
                self.assertTrue(done.wait(5))
                self.assertEqual(shared.ram_reserved, 60 * GIB)
                a.release_ram()
                b.release_ram()
                self.assertEqual(shared.ram_reserved, 0)


class RamWorkTests(unittest.TestCase):
    def test_ram_run_keeps_only_bam_and_qc_on_disk(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ram = root / 'ram'
            ram.mkdir()
            engine = Pipeline(root / 'data', root / 'refs', 1, lambda **kw: None, root / 'log.txt', ram=ram)
            engine.ram_root = ram / 'job'
            index = root / 'index' / 'genome'
            atomic_json(index.parent / 'complete.json', {'reference': 'mm10'})
            row = {'run_accession': 'SRR1', 'library_layout': 'SINGLE'}
            seen = {}

            def fake_map(row, index, directory, bam, signature, versions):
                # Inside the RAM work folder: SRA/FASTQ/fastp live there, the BAM lands in the data folder.
                seen['work'] = engine.work
                (engine.work / 'sra').mkdir(parents=True)
                (engine.work / 'sra' / 'SRR1.sra').write_bytes(b'sra')
                engine.ram_held = 5
                engine.shared.ram_reserved = 5
                bam.write_bytes(b'bam')
                return bam
            with patch.object(engine, 'map_reads', side_effect=fake_map):
                bam = engine.mapped_run(row, index, {'hisat2': 'h', 'samtools': 's', 'fastp': 'f', 'fasterq-dump': 'q'})
            self.assertEqual(seen['work'], ram / 'job' / 'SRR1')
            self.assertTrue(bam.is_relative_to(root / 'data' / 'bam'))
            self.assertFalse((ram / 'job' / 'SRR1').exists())
            self.assertFalse((root / 'data' / 'sra').exists())
            self.assertEqual(engine.shared.ram_reserved, 0)

    def test_ram_request_without_tmpfs_fails_clearly(self):
        with tempfile.TemporaryDirectory() as temp:
            engine = Pipeline(temp, temp, 1, lambda **kw: None, Path(temp) / 'log.txt', ram=None)
            with self.assertRaisesRegex(RuntimeError, 'RAM作業領域'):
                engine.run({'gene': 'X', 'assembly': 'mm10', 'wt': ['SRX1'], 'ko': ['SRX2'], 'ram_work': True}, Path(temp))


class ParallelBatchTests(unittest.TestCase):
    def test_batch_runs_up_to_parallel_items_at_once(self):
        running, peak, lock = [], [0], threading.Lock()
        release = threading.Event()

        def slow(engine, request, directory):
            with lock:
                running.append(request['gene'])
                peak[0] = max(peak[0], len(running))
            release.wait(5)
            with lock:
                running.remove(request['gene'])
            completed(engine, request, directory)
        text = 'gene,wt,ko\n' + ''.join(f'G{i},SRX{i}1,SRX{i}2\n' for i in range(5))
        with tempfile.TemporaryDirectory() as temp, patch.object(Pipeline, 'run', slow):
            jobs = Jobs(Path(temp) / 'jobs', Path(temp) / 'data', Path(temp) / 'refs', 1, parallel=3)
            record = jobs.create_batch({'csv': text, 'selected': list(range(5)), 'min_free_gb': 0, 'parallel': 2})
            time.sleep(0.5)
            view = jobs.batch_view(record['id'])
            self.assertEqual(view['counts'], {'running': 2, 'pending': 3})
            release.set()
            jobs.queue.join()
            self.assertEqual(jobs.batch_view(record['id'])['counts'], {'completed': 5})
            self.assertEqual(peak[0], 2)
            for value in (0, 4, 'x', True):
                with self.subTest(value=value), self.assertRaisesRegex(ValueError, '同時に実行'):
                    jobs.create_batch({'csv': text, 'selected': [0], 'parallel': value})

    def test_pausing_a_parallel_batch_pauses_all_its_running_analyses(self):
        started = threading.Semaphore(0)
        release = threading.Event()

        def slow(engine, request, directory):
            started.release()
            while not release.is_set():
                engine.checkpoint()
                time.sleep(0.05)
            completed(engine, request, directory)
        text = 'gene,wt,ko\n' + ''.join(f'G{i},SRX{i}1,SRX{i}2\n' for i in range(4))
        with tempfile.TemporaryDirectory() as temp, patch.object(Pipeline, 'run', slow):
            jobs = Jobs(Path(temp) / 'jobs', Path(temp) / 'data', Path(temp) / 'refs', 1, parallel=2)
            record = jobs.create_batch({'csv': text, 'selected': list(range(4)), 'min_free_gb': 0})
            self.assertTrue(started.acquire(timeout=5) and started.acquire(timeout=5))
            jobs.batch_control(record['id'], 'pause')
            jobs.queue.join()
            self.assertEqual(jobs.batch_view(record['id'])['counts'], {'paused': 2, 'pending': 2})
            release.set()
            jobs.batch_control(record['id'], 'resume')
            jobs.queue.join()
            self.assertEqual(jobs.batch_view(record['id'])['counts'], {'completed': 4})


if __name__ == '__main__':
    unittest.main()
