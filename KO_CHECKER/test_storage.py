"""Destructive-operation regressions; every file is disposable temporary data."""
import hashlib
import json
from pathlib import Path
import shutil
import stat
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

from storage import Storage, StorageBusy, is_link
from webapp import Jobs, handler


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.store = Storage(self.root / 'workdata', self.root / 'references')

    def file(self, path, data=b'test'):
        path = self.root / path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def entry(self, kind):
        return next(e for e in self.store.inventory()['entries'] if e['kind'] == kind)

    def test_inventory_sizes_only_generated_caches_and_does_not_delete(self):
        files = [self.file(p) for p in ('workdata/sra/SRR1/SRR1.sra', 'workdata/fastq/SRR1/read.fastq',
            'workdata/trimmed/SRR1/version/read.gz', 'workdata/bam/SRR1/version/read.bam',
            'workdata/samples/SRX1/version/sample.bam', 'references/mm10/genome.fa',
            'references/hisat2/' + 'a' * 24 + '/genome.1.ht2')]
        self.file('workdata/metadata/SRX1.json')
        self.file('workdata/unrelated/notes.txt')
        self.file('results/web/job/output/samples.csv')
        result = self.store.inventory()
        self.assertEqual(len(result['entries']), 7)
        self.assertEqual(result['total_bytes'], 28)
        self.assertTrue(all(p.exists() for p in files))
        self.assertTrue(all('_path' not in e and '_signature' not in e for e in result['entries']))

    def test_preview_and_delete_remove_companions_but_preserve_results(self):
        bam = self.file('workdata/bam/SRR1/version/read.bam')
        self.file('workdata/bam/SRR1/version/read.bam.bai')
        self.file('workdata/bam/SRR1/version/complete.json')
        result = self.file('results/web/job/output/samples.csv')
        raw = self.file('workdata/fastq/SRR1/read.fastq')
        preview = self.store.preview([self.entry('bam')['id']])
        self.assertTrue(bam.exists())
        self.assertEqual(preview['total_bytes'], 12)
        deleted = self.store.delete(preview['token'], True)
        self.assertEqual(deleted['deleted_bytes'], 12)
        self.assertFalse(bam.parent.parent.exists())
        self.assertTrue(result.exists() and raw.exists())
        self.assertEqual(deleted['failed'], [])

    def test_reference_and_index_can_be_deleted_without_removing_roots(self):
        self.file('references/mm10/genome.fa')
        self.file('references/mm10/genome.fa.fai')
        self.file('references/hisat2/' + 'a' * 24 + '/genome.1.ht2')
        preview = self.store.preview([e['id'] for e in self.store.inventory()['entries']])
        self.store.delete(preview['token'], True)
        self.assertTrue((self.root / 'references').is_dir())
        self.assertFalse((self.root / 'references/mm10').exists())
        self.assertEqual(self.store.inventory()['total_bytes'], 0)

    def test_confirmation_is_explicit_and_token_is_one_use(self):
        bam = self.file('workdata/bam/SRR1/read.bam')
        token = self.store.preview([self.entry('bam')['id']])['token']
        for confirmation in (False, None, 'true', 1):
            with self.assertRaises(ValueError):
                self.store.delete(token, confirmation)
            self.assertTrue(bam.exists())
        self.store.delete(token, True)
        with self.assertRaises(ValueError):
            self.store.delete(token, True)

    def test_expired_preview_cannot_delete(self):
        bam = self.file('workdata/bam/SRR1/read.bam')
        token = self.store.preview([self.entry('bam')['id']])['token']
        self.store.previews[token]['expires'] = -1
        with self.assertRaisesRegex(ValueError, '期限切れ'):
            self.store.delete(token, True)
        self.assertTrue(bam.exists())

    def test_changes_in_any_target_abort_before_any_deletion(self):
        bam = self.file('workdata/bam/SRR1/read.bam')
        fastq = self.file('workdata/fastq/SRR1/read.fastq')
        preview = self.store.preview([self.entry('bam')['id'], self.entry('fastq')['id']])
        fastq.write_bytes(b'changed after confirmation')
        with self.assertRaisesRegex(ValueError, '変更'):
            self.store.delete(preview['token'], True)
        self.assertTrue(bam.exists() and fastq.exists())

    def test_paths_unknown_ids_duplicates_and_invalid_inputs_are_rejected(self):
        bam = self.file('workdata/bam/SRR1/read.bam')
        key = self.entry('bam')['id']
        for ids in (None, [], True, [1], ['../results'], ['f' * 24], [key, key], [key] * 301):
            with self.subTest(ids=str(ids)[:40]), self.assertRaises(ValueError):
                self.store.preview(ids)
        self.assertTrue(bam.exists())

    def test_nested_symlink_is_never_followed_or_deleted(self):
        outside = self.file('outside/keep.txt', b'keep this')
        self.file('workdata/bam/SRR1/read.bam')
        (self.root / 'workdata/bam/SRR1/link').symlink_to(outside.parent, target_is_directory=True)
        entry = self.entry('bam')
        self.assertFalse(entry['deletable'])
        with self.assertRaises(ValueError):
            self.store.preview([entry['id']])
        self.assertEqual(outside.read_bytes(), b'keep this')

    def test_symlink_ancestor_and_windows_reparse_points_are_rejected(self):
        outside = self.file('outside/SRR1/keep.bam')
        (self.root / 'workdata').mkdir()
        (self.root / 'workdata/bam').symlink_to(outside.parent.parent, target_is_directory=True)
        result = self.store.inventory()
        self.assertTrue(result['warnings'])
        self.assertEqual(result['entries'], [])
        with self.assertRaises(ValueError):
            self.store.safe_path(self.store.data, self.store.data / 'bam/SRR1')
        self.assertTrue(outside.exists())
        self.assertTrue(is_link(SimpleNamespace(st_mode=stat.S_IFDIR, st_file_attributes=0x400)))

    def test_a_link_inserted_after_preview_cannot_escape(self):
        outside = self.file('outside/keep.txt')
        bam = self.file('workdata/bam/SRR1/read.bam')
        token = self.store.preview([self.entry('bam')['id']])['token']
        bam.unlink()
        bam.symlink_to(outside)
        with self.assertRaises(ValueError):
            self.store.delete(token, True)
        self.assertTrue(outside.exists())

    def test_delete_errors_are_reported_without_claiming_reclaimed_bytes(self):
        bam = self.file('workdata/bam/SRR1/read.bam')
        token = self.store.preview([self.entry('bam')['id']])['token']
        with patch('storage.shutil.rmtree', side_effect=PermissionError('in use')):
            result = self.store.delete(token, True)
        self.assertEqual(result['deleted_bytes'], 0)
        self.assertEqual(len(result['failed']), 1)
        self.assertTrue(bam.exists())

    def test_cleanup_blocks_submissions_until_files_are_removed(self):
        self.file('workdata/bam/SRR1/read.bam')
        with patch.object(Jobs, 'worker', lambda self: None):
            jobs = Jobs(self.root / 'results', self.store.data, self.store.references, 1)
        token = jobs.preview_storage([self.entry('bam')['id']])['token']
        entered, release, submitted = threading.Event(), threading.Event(), threading.Event()
        original = shutil.rmtree
        def removing(path):
            entered.set()
            release.wait(3)
            original(path)
        def submit():
            jobs.submit({'gene': 'Test', 'wt': 'SRX1', 'ko': 'SRX2'})
            submitted.set()
        with patch('storage.shutil.rmtree', side_effect=removing):
            deletion = threading.Thread(target=jobs.delete_storage, args=(token, True))
            deletion.start()
            self.assertTrue(entered.wait(2))
            submission = threading.Thread(target=submit)
            submission.start()
            try:
                self.assertFalse(submitted.wait(.05))
            finally:
                release.set()
                deletion.join(3)
                submission.join(3)
        self.assertTrue(submitted.is_set())


class StorageAPITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        with patch.object(Jobs, 'worker', lambda self: None):
            self.jobs = Jobs(self.root / 'results', self.root / 'data', self.root / 'refs', 1)
        self.bam = self.root / 'data/bam/SRR1/read.bam'
        self.bam.parent.mkdir(parents=True)
        self.bam.write_bytes(b'bam')
        self.report = self.root / 'results/saved.csv'
        self.report.write_bytes(b'results')
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), handler(self.jobs))
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.base = f'http://127.0.0.1:{self.server.server_port}'

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temp.cleanup()

    def get(self, path):
        with urllib.request.urlopen(self.base + path) as response:
            return json.load(response)

    def post(self, path, payload, headers=None):
        request = urllib.request.Request(self.base + path, data=json.dumps(payload).encode(),
            headers={'Content-Type': 'application/json', **(headers or {})})
        with urllib.request.urlopen(request) as response:
            return json.load(response)

    def test_http_preview_confirm_delete_and_origin_protection(self):
        inventory = self.get('/api/storage')
        self.assertFalse(inventory['busy'])
        key = inventory['entries'][0]['id']
        preview = self.post('/api/storage/preview', {'ids': [key]})
        self.assertTrue(self.bam.exists())
        for headers, status in (({'Origin': 'https://elsewhere.example'}, 403),
                                ({'Content-Type': 'text/plain'}, 415)):
            with self.assertRaises(urllib.error.HTTPError) as caught:
                self.post('/api/storage/delete', {'token': preview['token'], 'confirmed': True}, headers)
            self.assertEqual(caught.exception.code, status)
            self.assertTrue(self.bam.exists())
        with self.assertRaises(urllib.error.HTTPError) as caught:
            self.post('/api/storage/delete', {'token': preview['token'], 'confirmed': False})
        self.assertEqual(caught.exception.code, 400)
        result = self.post('/api/storage/delete', {'token': preview['token'], 'confirmed': True})
        self.assertEqual(result['deleted_bytes'], 3)
        self.assertFalse(self.bam.exists())
        self.assertEqual(self.report.read_bytes(), b'results')
        with urllib.request.urlopen(self.base + '/storage.js') as response:
            self.assertIn(b'/api/storage/preview', response.read())

    def test_running_and_queued_jobs_block_deletion_even_after_preview(self):
        key = self.get('/api/storage')['entries'][0]['id']
        for state in ('running', 'queued'):
            self.jobs.items['active'] = {'id': 'active', 'status': state}
            self.assertTrue(self.get('/api/storage')['busy'])
            with self.assertRaises(urllib.error.HTTPError) as caught:
                self.post('/api/storage/preview', {'ids': [key]})
            self.assertEqual(caught.exception.code, 409)
        self.jobs.items.clear()
        preview = self.post('/api/storage/preview', {'ids': [key]})
        self.jobs.items['active'] = {'id': 'active', 'status': 'running'}
        with self.assertRaises(urllib.error.HTTPError) as caught:
            self.post('/api/storage/delete', {'token': preview['token'], 'confirmed': True})
        self.assertEqual(caught.exception.code, 409)
        self.assertTrue(self.bam.exists())
        self.jobs.items['active']['status'] = 'paused'
        result = self.post('/api/storage/delete', {'token': preview['token'], 'confirmed': True})
        self.assertEqual(result['deleted_bytes'], 3)


if __name__ == '__main__':
    unittest.main()
