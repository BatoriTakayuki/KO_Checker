"""Local browser UI and persistent single-worker analysis queue (stdlib only)."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
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

from pipeline import Pipeline, PipelinePaused, ROOT, Shared, TOOLS, atomic_json, read_json, validate_request
from reporting import read_conclusion, conclusion_text
from storage import Storage, StorageBusy
import batch
import igv_view

FILES = {'samples.csv', 'events.csv', 'junctions.csv', 'consequences.csv', 'observations.csv', 'exon_coverage.csv', 'run.json',
         'region_deletions.csv', 'region_proteins.fasta', 'region_transcripts.fasta', 'gene_structure.json'}
# A batch stops dispatching after this many failures in a row (e.g. network outage, full disk).
BATCH_FAILURE_LIMIT = 5
SUMMARY_FIELDS = ['experiment_id', 'gene', 'status', 'sequence_status', 'ko_samples', 'ko_like', 'wt_like', 'undetermined',
                  'region_loss_ko_samples', 'wt_mean_cpm', 'ko_mean_cpm', 'ko_wt_ratio', 'sequence_text', 'overall',
                  'study', 'description', 'wt', 'ko', 'job_id', 'message']


def now():
    return datetime.now(timezone.utc).isoformat()


class Jobs:
    def __init__(self, root, data, references, threads, ram=None, parallel=1):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.data, self.references, self.threads = data, references, threads
        self.ram = ram
        # Analyses run at the same time; they share run/sample/reference locks and the RAM work area.
        self.parallel = max(1, int(parallel))
        self.shared = Shared(self.parallel)
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
        self.batch_root = self.root / '_batches'
        self.batch_root.mkdir(exist_ok=True)
        self.batches = {}
        for path in self.batch_root.glob('*.json'):
            record = read_json(path)
            if record:
                self.batches[record['id']] = record
        self.advance_batches()
        for _ in range(self.parallel):
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
            if sum(j['status'] in ('queued', 'running') and not j.get('batch_id') for j in self.items.values()) >= 10:
                raise ValueError('待機中の解析が多すぎます。完了後に再度実行してください。')
            return self.create(request)

    def create(self, request, **extra):
        with self.lock:
            jid = uuid.uuid4().hex
            job = {'id': jid, 'request': request, 'status': 'queued', 'stage': 'queued', 'progress': 0,
                   'threads': request.get('threads', self.threads),
                   'pause_requested': False,
                   'message': '実行待ち', 'created_at': now(), 'updated_at': now(), **extra}
            atomic_json(self.root / jid / 'job.json', job)
            self.items[jid] = job
            self.pause_events[jid] = threading.Event()
            self.queue.put(jid)
            return dict(job)

    # A batch keeps at most `parallel` of its items queued or running: the next item is queued only after
    # one has finished, so single analyses submitted meanwhile are not starved.
    def save_batch(self, record):
        record['updated_at'] = now()
        atomic_json(self.batch_root / (record['id'] + '.json'), record)

    def create_batch(self, payload):
        if not isinstance(payload, dict):
            raise ValueError('入力はJSONオブジェクトにしてください。')
        options = validate_request({'gene': 'X', 'wt': 'SRX1', 'ko': 'SRX2', 'assembly': payload.get('assembly', 'mm10'),
                                    'threads': payload.get('threads'), 'confirm_replicates': payload.get('confirm_replicates'),
                                    'ram_work': payload.get('ram_work')})
        multi_gene = payload.get('multi_gene', 'split')
        if multi_gene not in batch.MULTI_GENE:
            raise ValueError('複数遺伝子の行の扱いは split / first / skip から選択してください。')
        options.update(multi_gene=multi_gene)
        min_free = payload.get('min_free_gb', 50)
        if isinstance(min_free, str):
            min_free = min_free.strip() or 0
        try:
            min_free = float(min_free)
        except (TypeError, ValueError):
            min_free = -1
        if isinstance(payload.get('min_free_gb'), bool) or not 0 <= min_free <= 100000:
            raise ValueError('一時停止する空き容量は0〜100000 GBで指定してください。')
        parallel = payload.get('parallel', self.parallel)
        if isinstance(parallel, str):
            parallel = parallel.strip() or self.parallel
        if isinstance(parallel, bool) or not re.fullmatch(r'[0-9]+', str(parallel)) or not 1 <= int(parallel) <= self.parallel:
            raise ValueError(f'同時に実行する解析数は1〜{self.parallel}で指定してください（サーバーの同時実行数が上限です）。')
        rows = batch.parse_csv(payload.get('csv'))
        items = batch.build_items(rows, payload.get('selected'), options)
        name = str(payload.get('name') or 'batch.csv')[:120]
        with self.lock:
            if sum(b['status'] in ('running', 'paused') for b in self.batches.values()) >= 5:
                raise ValueError('未完了の一括処理が多すぎます。完了または取り消し後に作成してください。')
            bid = uuid.uuid4().hex
            record = {'id': bid, 'name': name, 'status': 'running', 'created_at': now(),
                      'options': {'assembly': options['assembly'], 'threads': options.get('threads'),
                                  'confirm_replicates': options.get('confirm_replicates') is True,
                                  'ram_work': options.get('ram_work') is True, 'parallel': int(parallel),
                                  'multi_gene': multi_gene, 'min_free_gb': min_free},
                      'items': items, 'failures_in_row': 0, 'message': '開始しました。'}
            self.batches[bid] = record
            self.save_batch(record)
            self.advance_batches()
            return self.batch_view(bid)

    def free_gb(self):
        path = Path(self.data) if Path(self.data).exists() else self.root
        return shutil.disk_usage(path).free / 1000**3

    def advance_batches(self):
        with self.lock:
            for record in self.batches.values():
                if record['status'] != 'running':
                    continue
                active = 0
                for item in record['items']:
                    if not item['job_id'] or item.get('settled'):
                        continue
                    job = self.items.get(item['job_id'])
                    status = job['status'] if job else 'failed'
                    if status in ('queued', 'running'):
                        active += 1
                    elif status == 'paused':
                        record.update(status='paused', message=f'{item["experiment_id"]} / {item["gene"]} の解析が中断されたため、一括処理を一時停止しました。「再開」で続きから実行します。')
                    else:
                        item['settled'] = True
                        record['failures_in_row'] = record['failures_in_row'] + 1 if status == 'failed' else 0
                if record['status'] == 'paused':
                    self.save_batch(record)
                    continue
                if record['failures_in_row'] >= BATCH_FAILURE_LIMIT:
                    record.update(status='paused', failures_in_row=0,
                                  message=f'{BATCH_FAILURE_LIMIT}件連続で失敗したため一時停止しました。ログで原因（通信・容量など）を確認してから「再開」してください。')
                    self.save_batch(record)
                    continue
                limit = min(record['options'].get('parallel', 1), self.parallel)
                minimum = record['options'].get('min_free_gb', 0)
                while active < limit:
                    pending = next((i for i, item in enumerate(record['items']) if not item['job_id'] and not item['cancelled']), None)
                    if pending is None:
                        if not active:
                            record.update(status='completed', message='すべての項目を処理しました。')
                        break
                    if minimum and self.free_gb() < minimum:
                        record.update(status='paused', message=f'保存先の空き容量が{minimum:g} GB未満のため一時停止しました。不要なデータを削除してから「再開」してください。')
                        break
                    item = record['items'][pending]
                    job = self.create(item['request'], batch_id=record['id'], batch_label=f'{item["experiment_id"]}')
                    item.update(job_id=job['id'], settled=False)
                    active += 1
                if record['status'] == 'running':
                    started = sum(1 for item in record['items'] if item['job_id'])
                    record['message'] = f'{started} / {len(record["items"])} 件目まで開始（同時実行 {active}件）。'
                self.save_batch(record)

    def batch_jobs(self, record, statuses):
        return [job for item in record['items'] if item['job_id'] and not item.get('settled')
                for job in [self.items.get(item['job_id'])] if job and job['status'] in statuses]

    def batch_control(self, bid, operation):
        with self.lock:
            record = self.batches[bid]
            if operation == 'pause':
                if record['status'] != 'running':
                    raise ValueError('実行中の一括処理だけを一時停止できます。')
                record.update(status='paused', message='利用者の操作で一時停止しました。')
                self.save_batch(record)
                for job in self.batch_jobs(record, ('queued', 'running')):
                    self.pause(job['id'])
            elif operation == 'resume':
                if record['status'] != 'paused':
                    raise ValueError('一時停止中の一括処理だけを再開できます。')
                record.update(status='running', message='再開しました。')
                self.save_batch(record)
                for job in self.batch_jobs(record, ('paused',)):
                    self.resume(job['id'])
            elif operation == 'cancel':
                if record['status'] not in ('running', 'paused'):
                    raise ValueError('未完了の一括処理だけを取り消せます。')
                for item in record['items']:
                    if not item['job_id']:
                        item['cancelled'] = True
                # Analyses already started are left as they are (they can be paused or resumed on their own).
                record.update(status='cancelled', message='未実行の項目を取り消しました。')
                self.save_batch(record)
            elif operation == 'retry':
                failed = [item for item in record['items']
                          if item['job_id'] and self.items.get(item['job_id'], {}).get('status') == 'failed']
                if not failed:
                    raise ValueError('失敗した項目はありません。')
                for item in failed:
                    item.update(job_id=None, cancelled=False, settled=False)
                if record['status'] == 'paused':
                    record['message'] = f'失敗した{len(failed)}件を未実行に戻しました。「再開」で実行します。'
                else:
                    record.update(status='running', failures_in_row=0, message=f'失敗した{len(failed)}件を再実行します。')
                self.save_batch(record)
            self.advance_batches()
            return self.batch_view(bid)

    def batch_item_status(self, item):
        if item['job_id']:
            job = self.items.get(item['job_id'])
            return job['status'] if job else 'failed'
        return 'cancelled' if item['cancelled'] else 'pending'

    def batch_view(self, bid, items=True):
        with self.lock:
            record = self.batches[bid]
            counts = {}
            rows = []
            for item in record['items']:
                status = self.batch_item_status(item)
                counts[status] = counts.get(status, 0) + 1
                if items:
                    job = self.items.get(item['job_id']) if item['job_id'] else None
                    rows.append({'experiment_id': item['experiment_id'], 'gene': item['gene'],
                                 'description': item.get('description', ''), 'fastq_gb': item.get('fastq_gb'),
                                 'wt': item['request']['wt'], 'ko': item['request']['ko'],
                                 'job_id': item['job_id'], 'status': status,
                                 'stage': job.get('stage', '') if job else '', 'progress': job.get('progress', 0) if job else 0,
                                 'message': job.get('message', '') if job else ''})
            view = {key: record[key] for key in ('id', 'name', 'status', 'created_at', 'updated_at', 'options', 'message')}
            view.update(total=len(record['items']), counts=counts)
            if items:
                view['items'] = rows
            return view

    def list_batches(self):
        with self.lock:
            return sorted((self.batch_view(bid, items=False) for bid in self.batches),
                          key=lambda b: b['created_at'], reverse=True)

    def batch_summary(self, bid):
        """One CSV row per batch item with the main conclusion of completed analyses."""
        with self.lock:
            record = self.batches[bid]
            entries = [(item, self.batch_item_status(item), dict(self.items.get(item['job_id']) or {})) for item in record['items']]
            confirmed = record['options'].get('confirm_replicates') is True
        handle = io.StringIO()
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for item, status, job in entries:
            row = {'experiment_id': item['experiment_id'], 'gene': item['gene'], 'status': status,
                   'study': item.get('study', ''), 'description': item.get('description', ''),
                   'wt': ' '.join(item['request']['wt']), 'ko': ' '.join(item['request']['ko']),
                   'job_id': item['job_id'] or '', 'message': job.get('message', '') if status != 'completed' else ''}
            if status == 'completed':
                try:
                    report = read_conclusion(self.root / item['job_id'] / 'output', confirmed)
                    counts = report['sequence']['counts']
                    expression = report['expression']
                    row.update(sequence_status=report['sequence']['status'], ko_samples=sum(counts.values()),
                               ko_like=counts['KO_like'], wt_like=counts['WT_like'], undetermined=counts['undetermined'],
                               region_loss_ko_samples=' '.join(report['region_loss']['ko_samples']),
                               wt_mean_cpm=expression.get('wt_mean_cpm'), ko_mean_cpm=expression.get('ko_mean_cpm'),
                               ko_wt_ratio=expression.get('ratio'), sequence_text=report['sequence']['text'],
                               overall=report['overall'])
                except (OSError, KeyError, ValueError) as error:
                    row['message'] = f'結論を読み取れません: {error}'
            writer.writerow({key: '' if value is None else value for key, value in row.items()})
        return ('﻿' + handle.getvalue()).encode('utf-8')

    def pause(self, jid):
        with self.lock:
            job = self.items[jid]
            if job['status'] == 'paused':
                return dict(job)
            if job['status'] not in ('queued', 'running'):
                raise ValueError('実行中・待機中の解析だけを中断できます。')
            if job['status'] == 'queued':
                self.update(jid, status='paused', message='実行前に中断しました。', paused_at=now())
                self.advance_batches()
            else:
                self.pause_events[jid].set()
                self.update(jid, pause_requested=True, message='中断処理中です。終了を確認するまでPCとDockerを起動したままにしてください。')
            return self.get(jid)

    def resume(self, jid):
        with self.lock:
            job = self.items[jid]
            if job['status'] != 'paused':
                raise ValueError('中断中の解析だけを再開できます。')
            if not job.get('batch_id') and sum(j['status'] in ('queued', 'running') and not j.get('batch_id')
                                               for j in self.items.values()) >= 10:
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
                                    lambda **kw: self.update(jid, **kw), directory / 'pipeline.log', event.is_set, self.ram, self.shared)
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
                try:
                    self.advance_batches()
                except Exception:
                    traceback.print_exc()
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
            batch_action = re.fullmatch(r'/api/batches(?:/(preview)|/([a-f0-9]{32})/(pause|resume|cancel|retry))?', self.path)
            if self.path != '/api/jobs' and not action and not storage_action and not batch_action:
                return self.send(404, {'error': 'Not found'})
            if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
                return self.send(415, {'error': 'application/jsonが必要です。'})
            try:
                length = int(self.headers.get('Content-Length', '0'))
                # Batch requests carry the uploaded CSV text.
                if not 0 < length <= (8 * 1024 * 1024 if batch_action else 16384):
                    raise ValueError('入力サイズが不正です。')
                payload = json.loads(self.rfile.read(length))
                if batch_action:
                    preview, bid, operation = batch_action.groups()
                    if preview:
                        if not isinstance(payload, dict):
                            raise ValueError('入力はJSONオブジェクトにしてください。')
                        rows = batch.parse_csv(payload.get('csv'))
                        return self.send(200, {'rows': rows, 'free_gb': round(jobs.free_gb(), 1)})
                    if bid:
                        if bid not in jobs.batches:
                            return self.send(404, {'error': '一括処理が見つかりません。'})
                        return self.send(200, jobs.batch_control(bid, operation))
                    return self.send(202, jobs.create_batch(payload))
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
            if path in ('/', '/app.js', '/storage.js', '/structure.js', '/batch.js', '/igv-panel.js', '/vendor/igv-3.0.0.min.js', '/style.css', '/stem-cell-boys.png'):
                file = ROOT / 'static' / ('index.html' if path == '/' else path[1:])
                return self.send(200, file.read_bytes(), mimetypes.guess_type(file.name)[0] + '; charset=utf-8')
            if path == '/api/health':
                names = TOOLS
                import rnako
                return self.send(200, {'tools': {n: bool(shutil.which(n)) for n in names}, 'pysam': rnako.pysam is not None,
                                       'default_threads': jobs.threads, 'logical_cpus': os.cpu_count(),
                                       'parallel': jobs.parallel,
                                       'ram_work_gb': (round(shutil.disk_usage(jobs.ram).total / 1024**3)
                                                       if jobs.ram and os.path.ismount(jobs.ram) else None),
                                       'disk_free_gb': round(shutil.disk_usage(jobs.root).free / 1024**3, 1)})
            if path == '/api/jobs':
                return self.send(200, jobs.list())
            if path == '/api/storage':
                return self.send(200, jobs.storage_inventory())
            if path == '/api/batches':
                return self.send(200, jobs.list_batches())
            batch_match = re.fullmatch(r'/api/batches/([a-f0-9]{32})(/summary\.csv)?', path)
            if batch_match:
                bid, summary = batch_match.groups()
                if bid not in jobs.batches:
                    return self.send(404, {'error': '一括処理が見つかりません。'})
                if not summary:
                    return self.send(200, jobs.batch_view(bid))
                content = jobs.batch_summary(bid)
                self.send_response(200)
                self.send_header('Content-Type', 'text/csv; charset=utf-8')
                self.send_header('Content-Disposition', f'attachment; filename="batch_{bid[:8]}_summary.csv"')
                self.send_header('Content-Length', str(len(content)))
                self.send_header('X-Content-Type-Options', 'nosniff')
                self.end_headers()
                self.wfile.write(content)
                return
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
    p.add_argument('--threads', type=int, default=int(os.environ.get('RNAKO_THREADS') or max(1, min(4, os.cpu_count() or 1))))
    p.add_argument('--parallel', type=int, default=int(os.environ.get('RNAKO_PARALLEL') or 1),
                   help='analyses run at the same time')
    p.add_argument('--ram', default='', help='tmpfs mount used when a job asks to work in RAM')
    args = p.parse_args()
    if args.threads < 1 or args.parallel < 1:
        p.error('--threads and --parallel must be positive')
    jobs = Jobs(args.jobs, args.data, args.references, args.threads, args.ram or None, args.parallel)
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
