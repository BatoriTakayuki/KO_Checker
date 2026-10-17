"""Real local HTTP and process tests for pause/resume; no public downloads."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from unittest.mock import patch
import urllib.request

from pipeline import Pipeline, PipelinePaused, validate_request
from transfers import download_file
from webapp import Jobs, handler


class DownloadServer:
    def __init__(self):
        self.data = b'ACGT' * 200000
        self.etag = '"v1"'
        self.ignore_range = False
        self.bad_range = False
        self.truncate_once = False
        self.validators = True
        self.requests = []
        owner = self
        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                owner.requests.append(dict(self.headers))
                requested = self.headers.get('Range')
                offset = int(requested.split('=')[1].split('-')[0]) if requested else 0
                use_range = bool(requested) and not owner.ignore_range and self.headers.get('If-Range') == owner.etag
                if use_range and offset >= len(owner.data):
                    self.send_response(416); self.end_headers(); return
                start = offset if use_range else 0
                body = owner.data[start:]
                self.send_response(206 if use_range else 200)
                self.send_header('Content-Length', str(len(body)))
                if owner.validators:
                    self.send_header('ETag', owner.etag)
                if use_range:
                    actual_start = 0 if owner.bad_range else start
                    self.send_header('Content-Range', f'bytes {actual_start}-{len(owner.data)-1}/{len(owner.data)}')
                self.end_headers()
                if owner.truncate_once:
                    owner.truncate_once = False
                    body = body[:len(body)//2]
                try:
                    self.wfile.write(body)
                except (BrokenPipeError, ConnectionResetError):
                    pass
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.url = f'http://127.0.0.1:{self.server.server_port}/reference.gz'

    def close(self):
        self.server.shutdown(); self.server.server_close(); self.thread.join()


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.server = DownloadServer()
        self.temp = tempfile.TemporaryDirectory()
        self.target = Path(self.temp.name) / 'reference.gz'

    def tearDown(self):
        self.server.close(); self.temp.cleanup()

    def pause_download(self):
        partial = self.target.with_name(self.target.name + '.partial')
        def check():
            if partial.exists() and partial.stat().st_size >= 128000:
                raise PipelinePaused('test pause')
        with self.assertRaises(PipelinePaused):
            download_file(self.server.url, self.target, check_pause=check)
        self.assertFalse(self.target.exists())
        self.assertGreater(partial.stat().st_size, 0)
        self.assertLess(partial.stat().st_size, len(self.server.data))
        return partial.stat().st_size

    def test_paused_download_resumes_using_range_and_validator(self):
        offset = self.pause_download()
        progress = []
        download_file(self.server.url, self.target, progress=lambda **kw: progress.append(kw))
        self.assertEqual(self.target.read_bytes(), self.server.data)
        self.assertEqual(self.server.requests[-1]['Range'], f'bytes={offset}-')
        self.assertEqual(self.server.requests[-1]['If-Range'], '"v1"')
        self.assertEqual(progress[0]['resumed_from'], offset)
        self.assertFalse(self.target.with_name(self.target.name + '.partial').exists())

    def test_server_ignoring_range_replaces_prefix_instead_of_appending(self):
        self.pause_download(); self.server.ignore_range = True
        download_file(self.server.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.server.data)

    def test_changed_remote_file_is_downloaded_fresh(self):
        self.pause_download()
        self.server.data = b'TGCA' * 150000
        self.server.etag = '"v2"'
        download_file(self.server.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.server.data)

    def test_incorrect_content_range_never_gets_appended(self):
        self.pause_download(); self.server.bad_range = True
        with patch('transfers.time.sleep'):
            download_file(self.server.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.server.data)
        self.assertNotIn('Range', self.server.requests[-1])

    def test_no_validator_means_safe_restart(self):
        self.server.validators = False
        self.pause_download()
        download_file(self.server.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.server.data)
        self.assertNotIn('Range', self.server.requests[-1])

    def test_truncated_transfer_retries_from_saved_prefix(self):
        self.server.truncate_once = True
        with patch('transfers.time.sleep'):
            download_file(self.server.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.server.data)
        self.assertIn('Range', self.server.requests[-1])

    def test_416_refetches_without_publishing_an_unverified_partial(self):
        partial = self.target.with_name(self.target.name + '.partial')
        partial.write_bytes(self.server.data)
        partial.with_name(partial.name + '.json').write_text(json.dumps({
            'url': self.server.url, 'validator': self.server.etag,
            'validator_header': 'ETag', 'total': len(self.server.data)}))
        with patch('transfers.time.sleep'):
            download_file(self.server.url, self.target)
        self.assertEqual(self.target.read_bytes(), self.server.data)
        self.assertNotIn('Range', self.server.requests[-1])


class ProcessTests(unittest.TestCase):
    def test_pause_terminates_running_command_and_can_run_next_command(self):
        with tempfile.TemporaryDirectory() as temp:
            event = threading.Event()
            pipeline = Pipeline(temp, temp, 1, lambda **kw: None, Path(temp)/'log.txt', event.is_set)
            timer = threading.Timer(0.3, event.set); timer.start()
            started = time.monotonic()
            try:
                with self.assertRaises(PipelinePaused):
                    pipeline.command([sys.executable, '-c', 'import time; time.sleep(60)'])
            finally:
                timer.cancel()
            self.assertLess(time.monotonic()-started, 8)
            event.clear()
            pipeline.command([sys.executable, '-c', 'print("resumed")'])
            self.assertIn('resumed', (Path(temp)/'log.txt').read_text())

    @unittest.skipUnless(os.name == 'posix', 'Linux process groups required')
    def test_pause_also_terminates_wrapper_child(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); child_pid = root/'child.pid'
            event = threading.Event()
            pipeline = Pipeline(temp,temp,1,lambda **kw:None,root/'log.txt',event.is_set)
            child_code = ('import signal,time,pathlib,os,sys; '
                          'signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                          'pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(60)')
            code = ('import subprocess,sys,time; '
                    f'p=subprocess.Popen([sys.executable,"-c",{child_code!r},sys.argv[1]]); '
                    'time.sleep(60)')
            timer = threading.Timer(0.5,event.set);timer.start()
            try:
                with self.assertRaises(PipelinePaused):
                    pipeline.command([sys.executable,'-c',code,str(child_pid)])
            finally:
                timer.cancel()
            self.assertTrue(child_pid.exists())
            pid = int(child_pid.read_text())
            # A killed child can briefly remain as a zombie awaiting init reaping.
            state = Path(f'/proc/{pid}/stat')
            for _ in range(20):
                if not state.exists() or state.read_text().split()[2] == 'Z':
                    break
                time.sleep(0.05)
            self.assertTrue(not state.exists() or state.read_text().split()[2] == 'Z')

    def test_reference_ready_pause_happens_before_any_run_is_downloaded(self):
        with tempfile.TemporaryDirectory() as temp:
            pipeline=Pipeline(temp,temp,1,lambda **kw:None,Path(temp)/'log.txt')
            row={'run_accession':'SRR1'}
            with patch.object(pipeline,'tools',return_value={}), patch('pipeline.resolve_runs',return_value=[row]), \
                 patch('pipeline.rnako.prepare_reference',return_value=(Path(temp)/'fa',Path(temp)/'gtf')), \
                 patch('pipeline.rnako.load_gene'),patch.object(pipeline,'index'),patch.object(pipeline,'mapped_run') as mapped:
                # Different experiments need different runs to pass duplicate checks.
                with patch('pipeline.resolve_runs',side_effect=[[row],[{'run_accession':'SRR2'}]]):
                    with self.assertRaisesRegex(PipelinePaused,'参照・索引'):
                        pipeline.run(validate_request({'gene':'G','wt':'SRX1','ko':'SRX2','pause_after_reference':True}),Path(temp))
                mapped.assert_not_called()


class JobControlTests(unittest.TestCase):
    def test_http_pause_resume_same_job_and_preserve_checkpoint(self):
        with tempfile.TemporaryDirectory() as temp:
            entered=threading.Event();calls=[]
            def run(engine,request,directory):
                calls.append(request)
                self.assertEqual(engine.threads,3)
                if len(calls)==1:
                    (directory/'checkpoint.txt').write_text('keep')
                    entered.set()
                    while True:
                        engine.checkpoint();time.sleep(0.02)
                self.assertEqual((directory/'checkpoint.txt').read_text(),'keep')
                engine.update(status='completed',stage='completed',progress=100,message='done')
            with patch.object(Pipeline,'run',run):
                jobs=Jobs(temp,temp,temp,1)
                server=ThreadingHTTPServer(('127.0.0.1',0),handler(jobs))
                thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
                try:
                    job=jobs.submit({'gene':'G','wt':'SRX1','ko':'SRX2','threads':3})
                    self.assertTrue(entered.wait(3))
                    def control(operation):
                        url=f'http://127.0.0.1:{server.server_port}/api/jobs/{job["id"]}/{operation}'
                        request=urllib.request.Request(url,data=b'{}',headers={'Content-Type':'application/json'})
                        with urllib.request.urlopen(request) as response:
                            self.assertEqual(response.status,202)
                            return json.load(response)
                    control('pause');jobs.queue.join()
                    self.assertEqual(jobs.get(job['id'])['status'],'paused')
                    jobs.threads=8  # Changing the server default must not change this job.
                    self.assertEqual(control('resume')['id'],job['id'])
                    jobs.queue.join()
                    self.assertEqual(jobs.get(job['id'])['status'],'completed')
                    self.assertEqual(len(calls),2)
                    self.assertEqual(jobs.get(job['id'])['threads'],3)
                finally:
                    server.shutdown();server.server_close();thread.join()

    def test_paused_queued_job_resumes_only_once(self):
        with tempfile.TemporaryDirectory() as temp:
            entered=threading.Event();release=threading.Event();calls=[]
            def run(engine,request,directory):
                calls.append(request['gene'])
                if request['gene']=='First':
                    entered.set();release.wait(3)
                engine.update(status='completed',stage='completed',message='done')
            with patch.object(Pipeline,'run',run):
                jobs=Jobs(temp,temp,temp,1)
                jobs.submit({'gene':'First','wt':'SRX1','ko':'SRX2'})
                self.assertTrue(entered.wait(3))
                job=jobs.submit({'gene':'Second','wt':'SRX3','ko':'SRX4'})
                self.assertEqual(jobs.pause(job['id'])['status'],'paused')
                jobs.resume(job['id']);release.set();jobs.queue.join()
                self.assertEqual(calls.count('Second'),1)

    def test_reference_ready_pause_does_not_repeat_after_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            calls=[]
            def run(engine,request,directory):
                calls.append(request)
                if request.get('pause_after_reference'):
                    engine.update(reference_ready=True)
                    raise PipelinePaused('ready')
                engine.update(status='completed',stage='completed',message='done')
            with patch.object(Pipeline,'run',run):
                jobs=Jobs(temp,temp,temp,1)
                job=jobs.submit({'gene':'G','wt':'SRX1','ko':'SRX2','pause_after_reference':True})
                jobs.queue.join();self.assertEqual(jobs.get(job['id'])['status'],'paused')
                jobs.resume(job['id']);jobs.queue.join()
                self.assertEqual(jobs.get(job['id'])['status'],'completed')
                self.assertNotIn('pause_after_reference',calls[1])


if __name__ == '__main__':
    unittest.main()
