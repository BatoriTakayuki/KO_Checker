"""Local browser UI and persistent single-worker analysis queue (stdlib only)."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
import os
from pathlib import Path
import queue
import re
import shutil
import threading
import traceback
import uuid
from urllib.parse import urlsplit

from pipeline import Pipeline, PipelinePaused, ROOT, atomic_json, read_json, validate_request

FILES = {'samples.csv', 'events.csv', 'junctions.csv', 'consequences.csv', 'observations.csv', 'exon_coverage.csv', 'run.json'}


def now():
    return datetime.now(timezone.utc).isoformat()


class Jobs:
    def __init__(self, root, data, references, threads):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.data, self.references, self.threads = data, references, threads
        self.lock = threading.RLock()
        self.queue = queue.Queue()
        self.items = {}
        self.pause_events = {}
        for path in self.root.glob('*/job.json'):
            job = read_json(path)
            if not job:
                continue
            if job['status'] in ('queued', 'running'):
                job.update(status='paused', stage=job.get('stage', 'interrupted'), pause_requested=False,
                           paused_at=now(), message='サーバーの再起動により中断しました。「再開」で完了済みデータを再利用します。', updated_at=now())
                atomic_json(path, job)
            self.items[job['id']] = job
        threading.Thread(target=self.worker, daemon=True).start()

    def list(self):
        with self.lock:
            return sorted((dict(j) for j in self.items.values()), key=lambda j: j['created_at'], reverse=True)

    def get(self, jid):
        with self.lock:
            job = self.items.get(jid)
            return dict(job) if job else None

    def update(self, jid, **changes):
        with self.lock:
            if 'stage' in changes and changes['stage'] != self.items[jid].get('stage'):
                self.items[jid].update(stage_started_at=now(), active_command='', activity='', transfer=None,
                                       storage_bytes=None, storage_label='')
            self.items[jid].update(**changes, updated_at=now())
            atomic_json(self.root / jid / 'job.json', self.items[jid])

    def submit(self, payload):
        request = validate_request(payload)
        with self.lock:
            if sum(j['status'] in ('queued', 'running') for j in self.items.values()) >= 10:
                raise ValueError('待機中の解析が多すぎます。完了後に再度実行してください。')
            jid = uuid.uuid4().hex
            job = {'id': jid, 'request': request, 'status': 'queued', 'stage': 'queued', 'progress': 0,
                   'threads': request.get('threads', self.threads),
                   'pause_requested': False,
                   'message': '実行待ち', 'created_at': now(), 'updated_at': now()}
            atomic_json(self.root / jid / 'job.json', job)
            self.items[jid] = job
            self.pause_events[jid] = threading.Event()
            self.queue.put(jid)
            return dict(job)

    def pause(self, jid):
        with self.lock:
            job = self.items[jid]
            if job['status'] == 'paused':
                return dict(job)
            if job['status'] not in ('queued', 'running'):
                raise ValueError('実行中・待機中の解析だけを中断できます。')
            if job['status'] == 'queued':
                self.update(jid, status='paused', message='実行前に中断しました。', paused_at=now())
            else:
                self.pause_events[jid].set()
                self.update(jid, pause_requested=True, message='中断処理中です。終了を確認するまでPCとDockerを起動したままにしてください。')
            return self.get(jid)

    def resume(self, jid):
        with self.lock:
            job = self.items[jid]
            if job['status'] != 'paused':
                raise ValueError('中断中の解析だけを再開できます。')
            if sum(j['status'] in ('queued', 'running') for j in self.items.values()) >= 10:
                raise ValueError('待機中の解析が多すぎます。')
            self.pause_events[jid] = threading.Event()
            self.update(jid, status='queued', pause_requested=False, message='再開待ち。完了済みデータを再利用します。')
            self.queue.put(jid)
            return self.get(jid)

    def worker(self):
        while True:
            jid = self.queue.get()
            try:
                with self.lock:
                    if self.items[jid]['status'] != 'queued':
                        continue
                    self.update(jid, status='running', attempt_started_at=now(), paused_at=None)
                    event = self.pause_events.setdefault(jid, threading.Event())
                    request = dict(self.items[jid]['request'])
                    threads = self.items[jid].get('threads', request.get('threads', self.threads))
                    self.update(jid, threads=threads)
                    if self.items[jid].get('reference_ready'):
                        request.pop('pause_after_reference', None)
                directory = self.root / jid
                pipeline = Pipeline(self.data, self.references, threads,
                                    lambda **kw: self.update(jid, **kw), directory / 'pipeline.log', event.is_set)
                pipeline.run(request, directory)
            except PipelinePaused as error:
                with (self.root / jid / 'pipeline.log').open('a', encoding='utf-8') as handle:
                    handle.write(str(error) + '\n')
                self.update(jid, status='paused', pause_requested=False, message=str(error), paused_at=now())
            except Exception as error:
                with (self.root / jid / 'pipeline.log').open('a', encoding='utf-8') as handle:
                    traceback.print_exc(file=handle)
                self.update(jid, status='failed', stage='failed', message=str(error))
            finally:
                self.queue.task_done()


def handler(jobs):
    class Handler(BaseHTTPRequestHandler):
        def send(self, status, data, content_type='application/json; charset=utf-8'):
            if not isinstance(data, bytes):
                data = json.dumps(data, ensure_ascii=False).encode('utf-8')
            self.send_response(status)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; object-src 'none'; frame-ancestors 'none'")
            self.end_headers()
            self.wfile.write(data)

        def do_POST(self):
            # Same-origin JSON requests only; reject cross-site form submissions.
            origin = self.headers.get('Origin')
            if (origin and urlsplit(origin).netloc != self.headers.get('Host')) or self.headers.get('Sec-Fetch-Site') == 'cross-site':
                return self.send(403, {'error': '同一オリジンから操作してください。'})
            action = re.fullmatch(r'/api/jobs/([a-f0-9]{32})/(pause|resume)', self.path)
            if self.path != '/api/jobs' and not action:
                return self.send(404, {'error': 'Not found'})
            if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
                return self.send(415, {'error': 'application/jsonが必要です。'})
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= 16384:
                    raise ValueError('入力サイズが不正です。')
                payload = json.loads(self.rfile.read(length))
                if action:
                    jid, operation = action.groups()
                    if not jobs.get(jid):
                        return self.send(404, {'error': '解析が見つかりません。'})
                    job = jobs.pause(jid) if operation == 'pause' else jobs.resume(jid)
                else:
                    job = jobs.submit(payload)
            except (ValueError, UnicodeError) as error:
                return self.send(400, {'error': str(error)})
            self.send(202, job)

        def do_GET(self):
            path = urlsplit(self.path).path
            if path in ('/', '/app.js', '/style.css'):
                file = ROOT / 'static' / ('index.html' if path == '/' else path[1:])
                return self.send(200, file.read_bytes(), mimetypes.guess_type(file.name)[0] + '; charset=utf-8')
            if path == '/api/health':
                names = ('prefetch', 'vdb-validate', 'fasterq-dump', 'hisat2', 'hisat2-build', 'samtools')
                import rnako
                return self.send(200, {'tools': {n: bool(shutil.which(n)) for n in names}, 'pysam': rnako.pysam is not None,
                                       'default_threads': jobs.threads, 'logical_cpus': os.cpu_count(),
                                       'disk_free_gb': round(shutil.disk_usage(jobs.root).free / 1024**3, 1)})
            if path == '/api/jobs':
                return self.send(200, jobs.list())
            match = re.fullmatch(r'/api/jobs/([a-f0-9]{32})(?:/(results|log|files/[^/]+))?', path)
            if not match or not jobs.get(match[1]):
                return self.send(404, {'error': '解析が見つかりません。'})
            jid, action = match.groups()
            job = jobs.get(jid)
            directory = jobs.root / jid
            if not action:
                return self.send(200, job)
            if action == 'log':
                log = directory / 'pipeline.log'
                if not log.exists():
                    return self.send(200, {'text': ''})
                with log.open('rb') as handle:
                    handle.seek(max(0, log.stat().st_size - 48000))
                    return self.send(200, {'text': handle.read().decode('utf-8', errors='replace')})
            if action == 'results':
                if job['status'] != 'completed':
                    return self.send(409, {'error': '解析はまだ完了していません。'})
                result = {'tables': {}, 'counts': {}}
                for name in ('samples', 'events', 'junctions', 'consequences', 'observations', 'exon_coverage'):
                    if name == 'exon_coverage' and not (directory / 'output' / (name + '.csv')).is_file():
                        result['tables'][name] = []
                        result['counts'][name] = 0
                        continue
                    # Bound browser payload; full CSV always available for download.
                    with (directory / 'output' / (name + '.csv')).open(newline='', encoding='utf-8') as handle:
                        rows = []
                        count = 0
                        for row in csv.DictReader(handle):
                            count += 1
                            if len(rows) < 1000:
                                rows.append(row)
                        result['tables'][name] = rows
                        result['counts'][name] = count
                return self.send(200, result)
            name = action.removeprefix('files/')
            if name in FILES and job['status'] == 'completed':
                file = directory / 'output' / name
            elif name in ('provenance.json', 'input_samples.csv', 'pipeline.log'):
                file = directory / name
            else:
                return self.send(404, {'error': 'ファイルが見つかりません。'})
            if not file.is_file():
                return self.send(404, {'error': 'ファイルはまだ作成されていません。'})
            self.send_response(200)
            self.send_header('Content-Type', 'application/octet-stream')
            self.send_header('Content-Disposition', f'attachment; filename="{name}"')
            self.send_header('Content-Length', str(file.stat().st_size))
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.end_headers()
            with file.open('rb') as handle:
                shutil.copyfileobj(handle, self.wfile)

    return Handler


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host', default='127.0.0.1')
    p.add_argument('--port', type=int, default=8080)
    p.add_argument('--data', default=str(ROOT / 'workdata'))
    p.add_argument('--references', default=str(ROOT / 'references'))
    p.add_argument('--jobs', default=str(ROOT / 'results' / 'web'))
    p.add_argument('--threads', type=int, default=max(1, min(4, os.cpu_count() or 1)))
    args = p.parse_args()
    if args.threads < 1:
        p.error('--threads must be positive')
    jobs = Jobs(args.jobs, args.data, args.references, args.threads)
    server = ThreadingHTTPServer((args.host, args.port), handler(jobs))
    print(f'KO_Checker Web: http://{args.host}:{args.port}', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == '__main__':
    main()
