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

from pipeline import Pipeline, PipelinePaused, ROOT, TOOLS, atomic_json, read_json, validate_request
from reporting import read_conclusion, conclusion_text
from storage import Storage, StorageBusy
import igv_view

FILES = {'samples.csv', 'events.csv', 'junctions.csv', 'consequences.csv', 'observations.csv', 'exon_coverage.csv', 'run.json',
         'region_deletions.csv', 'region_proteins.fasta', 'region_transcripts.fasta', 'gene_structure.json'}


def now():
    return datetime.now(timezone.utc).isoformat()


class Jobs:
    def __init__(self, root, data, references, threads):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.data, self.references, self.threads = data, references, threads
        self.lock = threading.RLock()
        self.storage = Storage(data, references)
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
            self.update(jid, status='queued', pause_requested=False, message='再開待ちです。完了済みデータを再利用します。')
            self.queue.put(jid)
            return self.get(jid)

    def storage_inventory(self):
        with self.lock:
            inventory = self.storage.inventory()
            active = [{'id': j['id'], 'gene': j.get('request', {}).get('gene', ''), 'status': j['status']}
                      for j in self.items.values() if j['status'] in ('running', 'queued')]
            return {**inventory, 'busy': bool(active), 'active_jobs': active}

    def storage_idle(self):
        if any(j['status'] in ('running', 'queued') for j in self.items.values()):
            raise StorageBusy('実行中・待機中の解析があるため削除できません。解析を完了するか、すべて中断してから操作してください。')

    def preview_storage(self, ids):
        with self.lock:
            self.storage_idle()
            return self.storage.preview(ids)

    def delete_storage(self, token, confirmed):
        with self.lock:
            self.storage_idle()
            return self.storage.delete(token, confirmed)

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
            self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; object-src 'none'; frame-ancestors 'none'")
            self.end_headers()
            self.wfile.write(data)

        def do_POST(self):
            # Same-origin JSON requests only; reject cross-site form submissions.
            origin = self.headers.get('Origin')
            if (origin and urlsplit(origin).netloc != self.headers.get('Host')) or self.headers.get('Sec-Fetch-Site') == 'cross-site':
                return self.send(403, {'error': '同一オリジンから操作してください。'})
            action = re.fullmatch(r'/api/jobs/([a-f0-9]{32})/(pause|resume)', self.path)
            storage_action = self.path in ('/api/storage/preview', '/api/storage/delete')
            if self.path != '/api/jobs' and not action and not storage_action:
                return self.send(404, {'error': 'Not found'})
            if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
                return self.send(415, {'error': 'application/jsonが必要です。'})
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= 16384:
                    raise ValueError('入力サイズが不正です。')
                payload = json.loads(self.rfile.read(length))
                if storage_action:
                    if not isinstance(payload, dict):
                        raise ValueError('入力はJSONオブジェクトにしてください。')
                    result = (jobs.preview_storage(payload.get('ids')) if self.path.endswith('/preview')
                              else jobs.delete_storage(payload.get('token'), payload.get('confirmed')))
                    return self.send(200, result)
                elif action:
                    jid, operation = action.groups()
                    if not jobs.get(jid):
                        return self.send(404, {'error': '解析が見つかりません。'})
                    job = jobs.pause(jid) if operation == 'pause' else jobs.resume(jid)
                else:
                    job = jobs.submit(payload)
            except StorageBusy as error:
                return self.send(409, {'error': str(error)})
            except (ValueError, UnicodeError, OSError) as error:
                return self.send(400, {'error': str(error)})
            self.send(202, job)

        def do_GET(self):
            path = urlsplit(self.path).path
            if path in ('/', '/app.js', '/storage.js', '/structure.js', '/igv-panel.js', '/vendor/igv-3.0.0.min.js', '/style.css', '/stem-cell-boys.png'):
                file = ROOT / 'static' / ('index.html' if path == '/' else path[1:])
                return self.send(200, file.read_bytes(), mimetypes.guess_type(file.name)[0] + '; charset=utf-8')
            if path == '/api/health':
                names = TOOLS
                import rnako
                return self.send(200, {'tools': {n: bool(shutil.which(n)) for n in names}, 'pysam': rnako.pysam is not None,
                                       'default_threads': jobs.threads, 'logical_cpus': os.cpu_count(),
                                       'disk_free_gb': round(shutil.disk_usage(jobs.root).free / 1024**3, 1)})
            if path == '/api/jobs':
                return self.send(200, jobs.list())
            if path == '/api/storage':
                return self.send(200, jobs.storage_inventory())
            igv_match = re.fullmatch(r'/api/jobs/([a-f0-9]{32})/igv/([^/]+)', path)
            if igv_match:
                jid, asset = igv_match.groups()
                job = jobs.get(jid)
                if not job or job['status'] != 'completed':
                    return self.send(404, {'error': 'Completed analysis not found'})
                output = jobs.root / jid / 'output'
                try:
                    if asset == 'config':
                        return self.send(200, igv_view.configuration(output, jid, jobs.data, jobs.references))
                    if asset == 'annotations.bed':
                        return self.send(200, igv_view.annotations(output)[1].encode(), 'text/plain; charset=utf-8')
                    files, _, _ = igv_view.assets(output, jobs.data, jobs.references)
                    if asset not in files:
                        return self.send(404, {'error': 'Unknown IGV asset'})
                    file = files[asset]
                    with file.open('rb') as handle:
                        size = os.fstat(handle.fileno()).st_size
                        start, end = 0, size-1
                        partial = self.headers.get('Range')
                        if partial:
                            try:
                                start, end = igv_view.byte_range(partial, size)
                            except ValueError:
                                self.send_response(416)
                                self.send_header('Content-Range', f'bytes */{size}')
                                self.send_header('Content-Length', '0')
                                self.end_headers()
                                return
                        self.send_response(206 if partial else 200)
                        self.send_header('Content-Type', 'application/octet-stream')
                        self.send_header('Accept-Ranges', 'bytes')
                        self.send_header('Content-Length', str(end-start+1))
                        self.send_header('Cache-Control', 'no-store')
                        if partial:
                            self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
                        self.end_headers()
                        handle.seek(start)
                        remaining = end-start+1
                        while remaining:
                            chunk = handle.read(min(1024*1024, remaining))
                            if not chunk:
                                break
                            self.wfile.write(chunk)
                            remaining -= len(chunk)
                    return
                except (ValueError, OSError, KeyError) as error:
                    return self.send(404, {'error': str(error)})
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
                result = {'tables': {}, 'counts': {}, 'region_analysis_available': (directory / 'output' / 'region_deletions.csv').is_file()}
                structure = directory / 'output' / 'gene_structure.json'
                result['structure'] = read_json(structure) if structure.is_file() else None
                for name in ('samples', 'events', 'junctions', 'consequences', 'observations', 'exon_coverage', 'region_deletions'):
                    if name in ('exon_coverage', 'region_deletions') and not (directory / 'output' / (name + '.csv')).is_file():
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
                result['conclusion'] = read_conclusion(directory / 'output', job.get('request', {}).get('confirm_replicates') is True)
                # Old CSVs used only WT_like baselines; present the current comparison
                # without altering previously saved analysis files.
                expression = result['conclusion']['expression']
                ratios = {s['srx_id']: s['ratio'] for s in expression['samples']}
                for row in result['tables']['samples']:
                    ratio = ratios.get(row.get('srx_id'))
                    row['expression_relative_to_wt'] = round(ratio, 5) if ratio is not None else ''
                    row['wt_baseline_cpm'] = expression['wt_mean_cpm'] if expression['wt_mean_cpm'] is not None else ''
                return self.send(200, result)
            name = action.removeprefix('files/')
            if name in ('conclusion.json', 'conclusion.txt') and job['status'] == 'completed':
                report = read_conclusion(directory / 'output', job.get('request', {}).get('confirm_replicates') is True)
                content = (json.dumps(report, ensure_ascii=False, indent=2) + '\n' if name.endswith('.json')
                           else conclusion_text(report)).encode('utf-8')
                self.send_response(200)
                self.send_header('Content-Type', 'application/octet-stream')
                self.send_header('Content-Disposition', f'attachment; filename="{name}"')
                self.send_header('Content-Length', str(len(content)))
                self.send_header('X-Content-Type-Options', 'nosniff')
                self.end_headers()
                self.wfile.write(content)
                return
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
