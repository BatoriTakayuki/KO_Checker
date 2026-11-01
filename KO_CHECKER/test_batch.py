"""Batch processing from an experiment list CSV."""
import csv
import io
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import batch
from pipeline import Pipeline
import webapp
from webapp import Jobs, handler

CSV = '''﻿WT・KOのファイル容量,,,,
"GB = 1,000,000,000 bytes",,,,
,,,,
experiment_id,target_label,total_FASTQ_GB,status,WT_SRX,KO_SRX,study,description,target_symbols
8162,Fus,0.046,取得済み,SRX1 SRX2,SRX3 SRX4,PRJNA1,Fus KO,Fus
3440,Gsk3a / Gsk3b,"1,112.9",取得済み,SRX5,SRX6,PRJNA2,DKO,Gsk3a / Gsk3b
603,Mettl3,,FASTQ未公開,SRX7,SRR8,PRJEB3,run only,Mettl3
'''


def completed(engine, request, directory):
    output = directory / 'output'
    output.mkdir()
    (output / 'samples.csv').write_text('srx_id,metadata_label,predicted_class,gene_unique_fragment_cpm,unique_primary_library_fragments\n'
                                        'SRX1,WT,undetermined,10,1000\nSRX3,KO,KO_like,4,1000\n', encoding='utf-8')
    engine.update(status='completed', stage='completed', progress=100, message='done')


class ParseTests(unittest.TestCase):
    def test_finds_header_below_title_lines_and_splits_genes(self):
        rows = batch.parse_csv(CSV)
        self.assertEqual([r['experiment_id'] for r in rows], ['8162', '3440', '603'])
        self.assertEqual(rows[0]['wt'], 'SRX1 SRX2')
        self.assertEqual(rows[1]['genes'], ['Gsk3a', 'Gsk3b'])
        self.assertEqual(rows[1]['fastq_gb'], 1112.9)
        self.assertIsNone(rows[2]['fastq_gb'])
        self.assertEqual(rows[0]['error'], '')
        self.assertIn('run番号', rows[2]['error'])

    def test_minimal_columns(self):
        rows = batch.parse_csv('gene,wt,ko\nPhf3,SRX1,SRX2\n')
        self.assertEqual(rows[0]['experiment_id'], '行2')
        self.assertEqual(rows[0]['genes'], ['Phf3'])

    def test_rejects_missing_columns_and_empty_input(self):
        for text in ('', 'a,b\n1,2\n', 'wt_srx,ko_srx\nSRX1,SRX2\n', 'gene,wt,ko\n'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                batch.parse_csv(text)

    def test_builds_one_item_per_gene_and_refuses_invalid_rows(self):
        rows = batch.parse_csv(CSV)
        items = batch.build_items(rows, [0, 1], {'assembly': 'mm39', 'threads': '8'})
        self.assertEqual([(i['experiment_id'], i['gene']) for i in items], [('8162', 'Fus'), ('3440', 'Gsk3a'), ('3440', 'Gsk3b')])
        self.assertEqual(items[0]['request'], {'gene': 'Fus', 'assembly': 'mm39', 'wt': ['SRX1', 'SRX2'],
                                               'ko': ['SRX3', 'SRX4'], 'threads': 8})
        first_only = batch.build_items(rows, [1], {'multi_gene': 'first'})
        self.assertEqual([i['gene'] for i in first_only], ['Gsk3a'])
        skipped = batch.build_items(rows, [0, 1], {'multi_gene': 'skip'})
        self.assertEqual([i['gene'] for i in skipped], ['Fus'])
        with self.assertRaisesRegex(ValueError, 'スキップ'):
            batch.build_items(rows, [1], {'multi_gene': 'skip'})
        for selected in ([2], [], [99], ['0'], [True]):
            with self.subTest(selected=selected), self.assertRaises(ValueError):
                batch.build_items(rows, selected, {})


class BatchJobTests(unittest.TestCase):
    def jobs(self, temp):
        return Jobs(Path(temp) / 'jobs', Path(temp) / 'data', Path(temp) / 'refs', 1)

    def test_runs_items_in_order_continues_after_failure_and_summarizes(self):
        calls = []

        def run(engine, request, directory):
            calls.append(request['gene'])
            if request['gene'] == 'Gsk3a':
                raise RuntimeError('gene not found')
            completed(engine, request, directory)
        with tempfile.TemporaryDirectory() as temp, patch.object(Pipeline, 'run', run):
            jobs = self.jobs(temp)
            record = jobs.create_batch({'csv': CSV, 'selected': [0, 1], 'assembly': 'mm39', 'min_free_gb': 0})
            jobs.queue.join()
            view = jobs.batch_view(record['id'])
            self.assertEqual(calls, ['Fus', 'Gsk3a', 'Gsk3b'])
            self.assertEqual(view['status'], 'completed')
            self.assertEqual(view['counts'], {'completed': 2, 'failed': 1})
            self.assertTrue(all(jobs.get(i['job_id'])['batch_id'] == record['id'] for i in view['items']))
            summary = list(csv.DictReader(io.StringIO(jobs.batch_summary(record['id']).decode('utf-8-sig'))))
            self.assertEqual([r['status'] for r in summary], ['completed', 'failed', 'completed'])
            self.assertEqual(summary[0]['ko_like'], '1')
            self.assertEqual(summary[0]['ko_wt_ratio'], '0.4')
            self.assertIn('gene not found', summary[1]['message'])

            calls.clear()
            with patch.object(Pipeline, 'run', completed):
                jobs.batch_control(record['id'], 'retry')
                jobs.queue.join()
            self.assertEqual(jobs.batch_view(record['id'])['counts'], {'completed': 3})

    def test_consecutive_failures_and_low_disk_pause_the_batch(self):
        def fail(engine, request, directory):
            raise RuntimeError('ENA unavailable')
        text = 'gene,wt,ko\n' + ''.join(f'G{i},SRX{i}1,SRX{i}2\n' for i in range(8))
        with tempfile.TemporaryDirectory() as temp, patch.object(Pipeline, 'run', fail):
            jobs = self.jobs(temp)
            record = jobs.create_batch({'csv': text, 'selected': list(range(8)), 'min_free_gb': 0})
            jobs.queue.join()
            view = jobs.batch_view(record['id'])
            self.assertEqual(view['status'], 'paused')
            self.assertEqual(view['counts'], {'failed': webapp.BATCH_FAILURE_LIMIT, 'pending': 8 - webapp.BATCH_FAILURE_LIMIT})
            with patch.object(Jobs, 'free_gb', return_value=10):
                low = jobs.create_batch({'csv': text, 'selected': [0], 'min_free_gb': 50})
            self.assertEqual(jobs.batch_view(low['id'])['status'], 'paused')
            self.assertIn('空き容量', jobs.batch_view(low['id'])['message'])

    def test_pause_resume_cancel_and_restart(self):
        started, release = threading.Event(), threading.Event()

        def slow(engine, request, directory):
            started.set()
            release.wait(5)
            engine.checkpoint()
            completed(engine, request, directory)
        with tempfile.TemporaryDirectory() as temp, patch.object(Pipeline, 'run', slow):
            jobs = self.jobs(temp)
            record = jobs.create_batch({'csv': CSV, 'selected': [0, 1], 'min_free_gb': 0})
            self.assertTrue(started.wait(5))
            jobs.batch_control(record['id'], 'pause')
            release.set()
            jobs.queue.join()
            view = jobs.batch_view(record['id'])
            self.assertEqual(view['status'], 'paused')
            self.assertEqual(view['counts'], {'paused': 1, 'pending': 2})
            with self.assertRaises(ValueError):
                jobs.batch_control(record['id'], 'pause')

            # A restarted server keeps the batch paused until the user resumes it.
            restarted = self.jobs(temp)
            self.assertEqual(restarted.batch_view(record['id'])['status'], 'paused')
            restarted.batch_control(record['id'], 'resume')
            restarted.queue.join()
            self.assertEqual(restarted.batch_view(record['id'])['counts'], {'completed': 3})

            started.clear()
            second = restarted.create_batch({'csv': CSV, 'selected': [0, 1], 'min_free_gb': 0})
            restarted.queue.join()
            self.assertEqual(restarted.batch_view(second['id'])['status'], 'completed')
            third = restarted.create_batch({'csv': CSV, 'selected': [0, 1], 'min_free_gb': 0})
            restarted.batch_control(third['id'], 'cancel')
            restarted.queue.join()
            view = restarted.batch_view(third['id'])
            self.assertEqual(view['status'], 'cancelled')
            self.assertEqual(view['counts'], {'completed': 1, 'cancelled': 2})

    def test_batch_jobs_do_not_count_against_single_job_limit(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(Jobs, 'worker'):
            jobs = self.jobs(temp)
            for _ in range(10):
                jobs.create({'gene': 'X', 'assembly': 'mm10', 'wt': ['SRX1'], 'ko': ['SRX2']}, batch_id='b' * 32)
            jobs.submit({'gene': 'Phf3', 'wt': 'SRX1', 'ko': 'SRX2'})


class BatchHttpTests(unittest.TestCase):
    def test_preview_create_list_and_summary_over_http(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(Pipeline, 'run', completed):
            jobs = Jobs(Path(temp) / 'jobs', Path(temp) / 'data', Path(temp) / 'refs', 1)
            server = ThreadingHTTPServer(('127.0.0.1', 0), handler(jobs))
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            base = f'http://127.0.0.1:{server.server_port}'

            def post(path, body):
                request = urllib.request.Request(base + path, data=json.dumps(body).encode(),
                                                 headers={'Content-Type': 'application/json'})
                with urllib.request.urlopen(request) as response:
                    return response.status, json.load(response)
            try:
                status, preview = post('/api/batches/preview', {'csv': CSV + 'x,y\n' * 3000})
                self.assertEqual(status, 200)
                self.assertEqual(preview['rows'][1]['genes'], ['Gsk3a', 'Gsk3b'])
                status, record = post('/api/batches', {'csv': CSV, 'selected': [0], 'name': 'list.csv', 'min_free_gb': '0'})
                self.assertEqual(status, 202)
                jobs.queue.join()
                with urllib.request.urlopen(base + '/api/batches') as response:
                    listed = json.load(response)
                    self.assertEqual(listed[0]['counts'], {'completed': 1})
                    self.assertNotIn('items', listed[0])
                with urllib.request.urlopen(base + f'/api/batches/{record["id"]}') as response:
                    self.assertEqual(json.load(response)['items'][0]['gene'], 'Fus')
                with urllib.request.urlopen(base + f'/api/batches/{record["id"]}/summary.csv') as response:
                    self.assertIn('attachment', response.headers['Content-Disposition'])
                    self.assertIn('8162', response.read().decode('utf-8-sig'))
                with urllib.request.urlopen(base + '/batch.js') as response:
                    self.assertEqual(response.status, 200)
                status, skipped = post('/api/batches', {'csv': CSV, 'selected': [0, 1], 'multi_gene': 'skip', 'min_free_gb': 0})
                self.assertEqual((skipped['total'], skipped['options']['multi_gene']), (1, 'skip'))
                jobs.queue.join()
                for body in ({'csv': CSV, 'selected': [2]}, {'csv': CSV, 'selected': [0], 'multi_gene': 'all'}):
                    with self.assertRaises(urllib.error.HTTPError) as caught:
                        post('/api/batches', body)
                    self.assertEqual(caught.exception.code, 400)
                with self.assertRaises(urllib.error.HTTPError) as caught:
                    post('/api/batches', {'csv': CSV, 'selected': [2]})
                self.assertEqual(caught.exception.code, 400)
                with self.assertRaises(urllib.error.HTTPError) as caught:
                    post(f'/api/batches/{"c" * 32}/pause', {})
                self.assertEqual(caught.exception.code, 404)
            finally:
                server.shutdown()
                server.server_close()
                thread.join()


if __name__ == '__main__':
    unittest.main()
