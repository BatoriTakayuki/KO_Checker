"""SRX -> prefetch -> fasterq-dump -> fastp -> HISAT2 BAM -> RNA evidence."""
from __future__ import annotations

import contextlib
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import queue
import re
import signal
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

import rnako
from reporting import write_conclusion

PIPELINE_VERSION = 2
# Raw SRA/FASTQ conversion has not changed; retain version-1 download caches.
FASTQ_VERSION = 1
FASTP_VERSION = 1
MULTIMAP_MAX = 100
TOOLS = ('prefetch', 'vdb-validate', 'fasterq-dump', 'fastp', 'hisat2', 'hisat2-build', 'samtools')
ACCESSION = re.compile(r"(?:SRX|ERX|DRX)[0-9]+\Z")
RUN = re.compile(r"(?:SRR|ERR|DRR)[0-9]+\Z")
ROOT = Path(__file__).resolve().parent
GIB = 1024 ** 3
# RAM needed per run from FASTQ conversion until its BAM is saved: uncompressed FASTQ is about 9x a
# normal .sra and 12x an SRA Lite archive, plus fastp output and the sort's temporary files.
RAM_PER_ARCHIVE_BYTE = 15
RAM_MARGIN = 4 * GIB
# Once the raw FASTQ is deleted, mapping holds the trimmed .fastq.gz, two uncompressed read blocks and the
# sort's temporary files plus the BAM (with -k 100 about 1-1.5x the trimmed files).
RAM_PER_TRIMMED_BYTE = 4
RAM_MAPPING_MARGIN = 8 * GIB
# HISAT2 2.2.1 with -k 100 keeps growing while it runs (~10 GiB after 4M read pairs, >200 GiB for a whole
# 60-80M pair library in TE_Checker, which the kernel killed). It is started afresh for each block of reads.
# Within a run HISAT2 reuses splice sites found in earlier reads, so a few spliced reads may align differently.
MAPPING_CHUNK_READS = 4_000_000
HISAT2_MEMORY_FLOOR = 32 * GIB


def fastp_options(layout):
    return ['--detect_adapter_for_pe'] if layout == 'PAIRED' else []


def mapping_options(layout):
    options = ['-k', str(MULTIMAP_MAX), '--secondary']
    if layout == 'PAIRED':
        options += ['--no-mixed', '--no-discordant']
    return options


def available_memory():
    """Bytes the container may still use (cgroup limit or MemAvailable), or None if unknown."""
    values = []
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemAvailable:'):
                values.append(int(line.split()[1]) * 1024)
    except (OSError, ValueError, IndexError):
        pass
    try:
        limit = Path('/sys/fs/cgroup/memory.max').read_text().strip()
        used = int(Path('/sys/fs/cgroup/memory.current').read_text().strip())
        if limit != 'max':
            values.append(int(limit) - used)
    except (OSError, ValueError):
        pass
    return min(values) if values else None


def sort_resources(threads, memory=None):
    """samtools sort threads and per-thread memory (MiB) that fit beside HISAT2 (a quarter of what is available)."""
    threads = max(1, min(int(threads), 8))
    if memory is None:
        return threads, 768
    per_thread = int(memory / 4 / (threads + 1) / 1024**2)
    while threads > 1 and per_thread < 256:
        threads -= 1
        per_thread = int(memory / 4 / (threads + 1) / 1024**2)
    return threads, max(128, min(768, per_thread))


def hisat2_memory_limit(memory, workers):
    """Address-space cap for one HISAT2 process, so a runaway run fails alone instead of exhausting the server."""
    if memory is None or not shutil.which('prlimit'):
        return None
    return max(HISAT2_MEMORY_FLOOR, memory // max(1, workers))


class ReadChunks:
    """Splits gzip FASTQ (one file, or the two mates) into blocks of `size` reads for successive HISAT2 runs.

    The hisat2 wrapper accepts only regular files, so each block is written uncompressed. A background thread
    writes the next block while HISAT2 maps the current one (two sets of files take turns).
    """
    BLOCK_LINES = 40000

    def __init__(self, fastqs, folder, size=None):
        self.size = size or MAPPING_CHUNK_READS
        self.fastqs = fastqs
        self.sets = [[Path(folder) / f'chunk{s}_{n}.fq' for n in range(1, len(fastqs) + 1)] for s in (0, 1)]
        self.readers = []  # started by the first next(), so an unused instance holds no processes
        self.free, self.ready = queue.Queue(), queue.Queue()
        self.producer = None
        self.current = None
        self.reads = 0

    def produce(self):
        try:
            while True:
                number = self.free.get()
                if number is None:
                    return
                counts = []
                for reader, path in zip(self.readers, self.sets[number]):
                    remaining, lines = self.size * 4, 0
                    with path.open('wb') as handle:
                        while remaining:
                            block = list(itertools.islice(reader.stdout, min(self.BLOCK_LINES, remaining)))
                            if not block:
                                break
                            handle.writelines(block)
                            lines += len(block)
                            remaining -= len(block)
                    counts.append(lines)
                if len(set(counts)) != 1 or counts[0] % 4:
                    raise RuntimeError(f'FASTQの行数が不正です（{counts}行）。ペアの_1と_2が対応していません。')
                if not counts[0]:
                    for reader in self.readers:
                        if reader.wait():
                            raise RuntimeError(f'FASTQの展開に失敗しました（gzip終了コード {reader.returncode}）')
                    self.ready.put(None)
                    return
                self.ready.put((number, counts[0] // 4))
        except Exception as error:
            self.ready.put(error)

    def next(self):
        """Paths of the next block's FASTQ files (the previous block is done with), or None at the end."""
        if self.producer is None:
            self.readers = [subprocess.Popen(['gzip', '-dc', str(p)], stdout=subprocess.PIPE,
                                             stdin=subprocess.DEVNULL) for p in self.fastqs]
            self.free.put(0)
            self.free.put(1)
            self.producer = threading.Thread(target=self.produce, daemon=True)
            self.producer.start()
        elif self.current is not None:
            self.free.put(self.current)
        item = self.ready.get()
        if isinstance(item, Exception):
            raise item
        if item is None:
            self.current = None
            return None
        self.current, reads = item
        self.reads += reads
        return self.sets[self.current]

    def close(self):
        for reader in self.readers:
            if reader.poll() is None:
                reader.kill()
        if self.producer is not None:
            self.free.put(None)
            self.producer.join(timeout=60)
        for reader in self.readers:
            reader.wait()
            reader.stdout.close()
        for path in self.sets[0] + self.sets[1]:
            path.unlink(missing_ok=True)


class Shared:
    """State shared by the analyses one server runs at the same time.

    Analyses of different experiments often use the same SRX (a common WT) and always the same reference,
    so each run, sample and reference is prepared by one analysis at a time; the others wait and reuse it.
    """
    def __init__(self, workers=1):
        self.workers = max(1, int(workers))
        self._guard = threading.Lock()
        self._locks = {}
        self.ram = threading.Condition()
        self.ram_reserved = 0

    def lock(self, key):
        with self._guard:
            return self._locks.setdefault(key, threading.Lock())


class PipelinePaused(RuntimeError):
    pass


def stop_process(process):
    """Stop this command's process group, including HISAT2 wrapper children."""
    if os.name == 'posix':
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    elif process.poll() is None:
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        if os.name == 'posix':
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
    if os.name == 'posix':
        # A wrapper may exit on SIGTERM while a child ignores it. Do not leave
        # that child using resources after the job has been labelled paused.
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Unique per writer: analyses running at the same time may refresh the same cache file.
    temporary = path.with_suffix(path.suffix + f'.{os.getpid()}.{threading.get_ident()}.partial')
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
    os.replace(temporary, path)


def read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None


def fingerprint(paths):
    return [{'path': str(Path(p).resolve()), 'size': Path(p).stat().st_size,
             'mtime_ns': Path(p).stat().st_mtime_ns} for p in paths]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:24]


def valid_cache(directory, signature):
    manifest = read_json(Path(directory) / 'complete.json')
    if not manifest or manifest.get('signature') != signature:
        return False
    try:
        paths = [entry['path'] for entry in manifest['files']]
        return bool(paths) and fingerprint(paths) == manifest['files'] and all(Path(p).stat().st_size for p in paths)
    except (OSError, KeyError, TypeError):
        return False


def mark_complete(directory, signature, paths):
    atomic_json(Path(directory) / 'complete.json', {'signature': signature, 'files': fingerprint(paths)})


def validate_request(payload):
    if not isinstance(payload, dict):
        raise ValueError('入力はJSONオブジェクトにしてください。')
    gene = payload.get('gene', '')
    if not isinstance(gene, str) or not re.fullmatch(r'[A-Za-z0-9_.-]{1,120}', gene.strip()):
        raise ValueError('対象遺伝子は遺伝子名またはgene IDで指定してください。')
    assembly = payload.get('assembly', 'mm10')
    if assembly not in rnako.REFERENCE:
        raise ValueError('参照ゲノムはmm10 / mm39 / hg38から選択してください。')
    groups = {}
    for label in ('wt', 'ko'):
        value = payload.get(label, '')
        if not isinstance(value, str):
            raise ValueError('アクセッションは文字列で指定してください。')
        ids = [x.upper() for x in re.split(r'[\s,;、]+', value.strip()) if x]
        if not ids or len(ids) > 30 or any(not ACCESSION.fullmatch(x) for x in ids):
            raise ValueError('WT・KOそれぞれにSRX（ERX/DRXも可）を1〜30件指定してください。run番号は入力できません。')
        if len(ids) != len(set(ids)):
            raise ValueError('同じ群にアクセッションが重複しています。')
        groups[label] = ids
    if set(groups['wt']) & set(groups['ko']):
        raise ValueError('同じアクセッションをWTとKOの両方には指定できません。')
    result = {'gene': gene.strip(), 'assembly': assembly, **groups}
    threads = payload.get('threads')
    if isinstance(threads, str):
        threads = threads.strip()
    if threads is not None and threads != '':
        if (isinstance(threads, bool) or not isinstance(threads, (str, int))
                or not re.fullmatch(r'[0-9]+', str(threads)) or not 1 <= int(threads) <= 256):
            raise ValueError('CPUスレッド数は1〜256の整数で指定してください。空欄なら既定値を使用します。')
        result['threads'] = int(threads)
    if payload.get('pause_after_reference') is True:
        result['pause_after_reference'] = True
    if payload.get('confirm_replicates') is True:
        result['confirm_replicates'] = True
    if payload.get('ram_work') is True:
        result['ram_work'] = True
    return result


def fetch_run_metadata(accession, note=None):
    """Retry transient ENA failures, then use its independent search endpoint."""
    note = note or (lambda message: None)
    common = {'result': 'read_run',
              'fields': 'run_accession,experiment_accession,library_layout,library_strategy,tax_id,scientific_name',
              'format': 'json'}
    routes = [('filereport', {'accession': accession}),
              ('search', {'query': f'experiment_accession="{accession}"', 'limit': '0'})]
    errors = []
    successful_response = False
    for endpoint, params in routes:
        url = 'https://www.ebi.ac.uk/ena/portal/api/' + endpoint + '?' + urllib.parse.urlencode({**common, **params})
        for attempt in range(1, 4):
            note(f'{accession}: ENA {endpoint}でrun情報を取得中 ({attempt}/3)')
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'RNAKO-Web/1.0', 'Accept': 'application/json'})
                with urllib.request.urlopen(request, timeout=45) as response:
                    rows = json.load(response)
                if not isinstance(rows, list):
                    raise ValueError('応答のJSONがrun一覧ではありません。')
                successful_response = True
                if rows:
                    return rows
                note(f'{accession}: ENA {endpoint}では該当runが見つかりません。')
                break
            except urllib.error.HTTPError as error:
                description = f'HTTP {error.code}'
                retryable = error.code in (408, 429, 500, 502, 503, 504)
            except (urllib.error.URLError, TimeoutError, OSError, ValueError) as error:
                description = str(error)
                retryable = True
            errors.append(f'{endpoint}: {description}')
            note(f'{accession}: ENA {endpoint}の取得失敗: {description}')
            if not retryable or attempt == 3:
                break
            delay = 2 ** attempt
            note(f'{delay}秒後に再試行します。')
            time.sleep(delay)
        if endpoint == 'filereport':
            note(f'{accession}: ENA検索APIでの取得に切り替えます。')
    if successful_response and not errors:
        raise ValueError(f'{accession}: 公開runが見つかりません。番号や公開状態を確認してください。')
    raise RuntimeError(f'{accession}: 外部サービスENAからrun情報を取得できませんでした。'
                       '時間をおいて再実行してください。詳細: ' + '; '.join(errors[-3:]))


def resolve_runs(accession, assembly, cache, note=None):
    """ENA indexes the INSDC/SRA experiments, including their constituent runs."""
    if not ACCESSION.fullmatch(accession):
        raise ValueError('SRX（ERX/DRXも可）を指定してください。run番号は入力できません。')
    path = Path(cache) / 'metadata' / (accession + '.json')
    cached = read_json(path)
    if cached and time.time() - cached.get('fetched_at', 0) < 86400:
        rows = cached['runs']
        if note:
            note(f'再利用 SRXのrun情報: {accession}')
    else:
        rows = fetch_run_metadata(accession, note)
        atomic_json(path, {'fetched_at': time.time(), 'runs': rows})
    expected_taxon = '9606' if assembly == 'hg38' else '10090'
    if not isinstance(rows, list) or not rows:
        raise ValueError(f'{accession}: run情報の件数が不正です。')
    seen = set()
    for row in rows:
        run = row.get('run_accession', '')
        if not RUN.fullmatch(run) or run in seen or row.get('experiment_accession') != accession:
            raise ValueError(f'{accession}: run情報が不正です。')
        seen.add(run)
        if str(row.get('tax_id')) != expected_taxon:
            raise ValueError(f'{accession}: 生物種が{assembly}と一致しません ({row.get("scientific_name")})。')
        if row.get('library_strategy') != 'RNA-Seq' or row.get('library_layout') not in ('PAIRED', 'SINGLE'):
            raise ValueError(f'{accession}: 対応するRNA-Seqライブラリではありません。')
    return sorted(rows, key=lambda r: r['run_accession'])


class Pipeline:
    def __init__(self, data, references, threads, update, log, should_pause=None, ram=None, shared=None):
        self.data = Path(data).resolve()
        self.references = Path(references).resolve()
        self.threads = threads
        self._update = update
        self.log = Path(log)
        self.should_pause = should_pause or (lambda: False)
        self.last_heartbeat = 0
        self.current_run = None
        self.current_index = None
        self.current_trim = None
        self.reference_fasta = None
        self.preprocessing = {}
        self.ram = Path(ram) if ram else None
        self.ram_root = None
        # RAM work directory of the run being processed (None: work on disk).
        self.work_dir = None
        self.shared = shared or Shared()
        self.ram_held = 0

    @property
    def work(self):
        """Where SRA, FASTQ and fastp files are written: the run's RAM directory or the data folder."""
        return self.work_dir or self.data

    @contextlib.contextmanager
    def exclusive(self, key, waiting):
        """Hold a resource other analyses running at the same time may also prepare; wait (pausable) if busy."""
        lock = self.shared.lock(key)
        if not lock.acquire(blocking=False):
            self.note(waiting)
            self.update(activity=waiting)
            while not lock.acquire(timeout=1):
                self.checkpoint()
        try:
            yield
        finally:
            lock.release()

    def reserve_ram(self, run, amount):
        """Hold part of the RAM work area until this run's BAM is saved; wait while other analyses use it."""
        capacity = shutil.disk_usage(self.ram).total
        if amount > capacity * 0.9:
            raise RuntimeError(f'{run}: RAM作業領域（{capacity / GIB:.0f} GiB）に収まりません'
                               f'（見積もり {amount / GIB:.0f} GiB）。RAMを使わない設定で実行してください。')
        shared = self.shared
        with shared.ram:
            if shared.ram_reserved + amount > capacity * 0.9:
                self.note(f'{run}: RAM作業領域の空き待ち（見積もり {amount / GIB:.0f} GiB、'
                          f'使用予約 {shared.ram_reserved / GIB:.0f} / {capacity / GIB:.0f} GiB）')
            while shared.ram_reserved + amount > capacity * 0.9:
                shared.ram.wait(timeout=1)
                self.checkpoint()
            shared.ram_reserved += amount
        self.ram_held += amount

    def shrink_ram(self, amount):
        """Lower this run's reservation (never raise it) so waiting analyses can start converting."""
        if self.ram_held > amount:
            with self.shared.ram:
                self.shared.ram_reserved -= self.ram_held - amount
                self.shared.ram.notify_all()
            self.ram_held = amount

    def release_ram(self):
        self.shrink_ram(0)

    def checkpoint(self):
        if self.should_pause():
            raise PipelinePaused('利用者の操作で中断しました。完了済みデータを再利用して再開できます。')

    def update(self, **changes):
        self.checkpoint()
        self._update(**changes)

    def heartbeat(self, command):
        if time.monotonic() - self.last_heartbeat < 2:
            return
        activity = ''
        if self.log.exists():
            with self.log.open('rb') as handle:
                handle.seek(max(0, self.log.stat().st_size - 4096))
                lines = handle.read().decode('utf-8', errors='replace').splitlines()
                activity = next((line.strip()[:240] for line in reversed(lines) if line.strip()), '')
        storage_bytes = None
        storage_label = ''
        if command == 'prefetch' and self.current_run:
            folder = self.work / 'sra' / self.current_run
            paths = (p for p in folder.rglob('*') if p.is_file() and p.suffix not in ('.lock', '.prf'))
            storage_label = 'SRA保存領域'
        elif command == 'fasterq-dump' and self.current_run:
            paths = (self.work / 'fastq' / self.current_run).glob('*.fastq')
            storage_label = 'FASTQ書き出し済み'
        elif command == 'fastp' and self.current_trim:
            paths = self.current_trim.glob('*.partial.fastq.gz')
            storage_label = '前処理済みFASTQ書き出し済み'
        elif command == 'HISAT2 / samtools sort' and self.current_run:
            paths = ((self.work / 'bam' / self.current_run) if self.work_dir else (self.data / 'bam' / self.current_run)).rglob('*.partial.bam')
            storage_label = 'BAM書き出し済み（バッファ中は増えません）'
        elif command == 'hisat2-build' and self.current_index:
            paths = self.current_index.glob('genome.*.ht2*')
            storage_label = '索引ファイル書き出し済み'
        else:
            paths = None
        if paths is not None:
            try:
                storage_bytes = sum(p.stat().st_size for p in paths)
            except OSError:
                pass
        self.update(active_command=command, activity=activity,
                    storage_bytes=storage_bytes, storage_label=storage_label)
        self.last_heartbeat = time.monotonic()

    def wait_process(self, process, command):
        while process.poll() is None:
            self.checkpoint()
            self.heartbeat(command)
            time.sleep(0.25)
        return process.wait()

    def note(self, message):
        self.checkpoint()
        with self.log.open('a', encoding='utf-8') as handle:
            handle.write(message + '\n')

    def command(self, args, stdout=None):
        self.checkpoint()
        self.note('$ ' + ' '.join(map(str, args)))
        with self.log.open('ab') as log:
            process = subprocess.Popen(list(map(str, args)), stdout=stdout if stdout is not None else log,
                                       stderr=log, stdin=subprocess.DEVNULL, start_new_session=(os.name == 'posix'))
            try:
                returncode = self.wait_process(process, Path(str(args[0])).name)
            except BaseException:
                stop_process(process)
                raise
        if returncode:
            raise RuntimeError(f'{Path(str(args[0])).name} が終了コード{returncode}で失敗しました。ログを確認してください。')

    def tools(self):
        names = TOOLS
        missing = [name for name in names if not shutil.which(name)]
        if rnako.pysam is None or missing:
            raise RuntimeError('解析環境が不足しています。Dockerで起動してください。不足: ' + ', '.join(missing + ([] if rnako.pysam else ['pysam'])))
        versions = {}
        for name in names:
            result = subprocess.run([name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            versions[name] = result.stdout.decode('utf-8', errors='replace').strip()
        return versions

    def fastqs(self, row, versions):
        run = row['run_accession']
        self.current_run = run
        signature = {'version': FASTQ_VERSION, 'run': run, 'layout': row['library_layout'],
                     'fasterq': versions['fasterq-dump'], 'mode': 'split-files'}
        directory = self.work / 'fastq' / run
        directory.mkdir(parents=True, exist_ok=True)
        paired = row['library_layout'] == 'PAIRED'
        paths = [directory / (run + suffix + '.fastq') for suffix in (('_1', '_2') if paired else ('',))]
        if valid_cache(directory, signature):
            self.note(f'再利用 FASTQ: {run}')
            return paths
        archive_root = self.work / 'sra'
        archive_root.mkdir(parents=True, exist_ok=True)
        # prefetch resumes interrupted downloads and validates archive availability.
        self.update(stage='sra_download', message=f'{run}: SRAデータを取得中', transfer=None)
        self.command(['prefetch', run, '--max-size', 'u', '--output-directory', archive_root / run])
        # Depending on prefetch's output-directory handling, -O <root>/<run>
        # places the archive directly there or creates another <run> directory.
        # Only accept this accession's exact filename in these two known layouts.
        candidates = [parent / (run + suffix)
                      for parent in (archive_root / run, archive_root / run / run)
                      for suffix in ('.sra', '.sralite', '')]
        archive = next((p for p in candidates if p.is_file() and p.stat().st_size), None)
        if archive is None:
            raise RuntimeError(f'{run}: prefetch後のSRAファイルが見つかりません。確認先: '
                               + ', '.join(map(str, candidates)))
        self.update(stage='sra_validate', message=f'{run}: 取得データの整合性を検査中')
        self.command(['vdb-validate', archive])
        if self.work_dir:
            self.reserve_ram(run, archive.stat().st_size * RAM_PER_ARCHIVE_BYTE + RAM_MARGIN)
        temp = directory / 'scratch'
        temp.mkdir(exist_ok=True)
        self.update(stage='fastq_convert', message=f'{run}: SRAをFASTQへ変換中')
        self.command(['fasterq-dump', archive, '--split-files', '--force', '--threads', self.threads,
                      '--outdir', directory, '--temp', temp])
        if not paired and not paths[0].exists() and (directory / (run + '_1.fastq')).exists():
            os.replace(directory / (run + '_1.fastq'), paths[0])
        if not all(p.is_file() and p.stat().st_size for p in paths):
            raise RuntimeError(f'{run}: 期待するFASTQが揃いません。ライブラリ構成を確認してください。')
        unexpected = [p for p in directory.glob('*.fastq') if p not in paths and p.stat().st_size]
        if unexpected:
            raise RuntimeError(f'{run}: 追加のFASTQが存在します。特殊なライブラリ構成には未対応です。')
        mark_complete(directory, signature, paths)
        if self.work_dir:
            shutil.rmtree(archive_root / run)  # the RAM copy of the archive is not needed after conversion
            shutil.rmtree(temp, ignore_errors=True)
        return paths

    def trimmed_fastqs(self, row, raw, versions):
        run, layout = row['run_accession'], row['library_layout']
        paired = layout == 'PAIRED'
        if len(raw) != (2 if paired else 1):
            raise RuntimeError(f'{run}: fastpの入力FASTQ数がライブラリ構成と一致しません。')
        options = fastp_options(layout)
        signature = {'version': FASTP_VERSION, 'run': run, 'layout': layout,
                     'raw_fastq': fingerprint(raw), 'fastp': versions['fastp'], 'options': options}
        directory = self.work / 'trimmed' / run / digest(signature)
        directory.mkdir(parents=True, exist_ok=True)
        self.current_trim = directory
        paths = [directory / (run + suffix + '.trim.fastq.gz')
                 for suffix in (('_1', '_2') if paired else ('',))]
        html, report = directory / (run + '_fastp.html'), directory / (run + '_fastp.json')
        metadata = directory / 'preprocessing.json'
        # Debian fastp releases cap their worker count at 16.
        threads = min(self.threads, 16)
        if valid_cache(directory, signature):
            self.note(f'再利用 fastp前処理済みFASTQ: {run}')
        else:
            self.update(stage='fastq_trim', message=f'{run}: fastpでアダプタ除去・品質フィルタ中', transfer=None)
            partial = [p.with_name(p.name.replace('.trim.fastq.gz', '.trim.partial.fastq.gz')) for p in paths]
            partial_html, partial_report = directory / (run + '_fastp.partial.html'), directory / (run + '_fastp.partial.json')
            args = ['fastp', *options, '-i', raw[0], '-o', partial[0]]
            if paired:
                args += ['-I', raw[1], '-O', partial[1]]
            args += ['-h', partial_html, '-j', partial_report, '-w', threads]
            self.command(args)
            self.checkpoint()
            summary = read_json(partial_report)
            try:
                remaining = summary['summary']['after_filtering']['total_reads']
            except (TypeError, KeyError):
                raise RuntimeError(f'{run}: fastpのQCレポートが不完全です。') from None
            if not isinstance(remaining, int) or isinstance(remaining, bool) or remaining <= 0:
                raise RuntimeError(f'{run}: fastp処理後に有効なリードが残りません。QCレポートを確認してください。')
            sources = [*partial, partial_html, partial_report]
            if not all(p.is_file() and p.stat().st_size for p in sources):
                raise RuntimeError(f'{run}: fastpの出力FASTQ・QCレポートが揃いません。')
            for source, destination in zip(sources, [*paths, html, report]):
                os.replace(source, destination)
            atomic_json(metadata, {'raw_fastq': list(map(str, raw)), 'trimmed_fastq': list(map(str, paths)),
                                   'fastp_html': str(html), 'fastp_json': str(report),
                                   'fastp_options': options, 'fastp_threads': threads})
            mark_complete(directory, signature, [*paths, html, report, metadata])
        self.preprocessing[run] = read_json(metadata)
        return paths

    def index(self, fasta, versions):
        self.reference_fasta = Path(fasta)
        signature = {'fasta': fingerprint([fasta]), 'hisat2-build': versions['hisat2-build']}
        directory = self.references / 'hisat2' / digest(signature)
        directory.mkdir(parents=True, exist_ok=True)
        prefix = directory / 'genome'
        self.current_index = directory
        if valid_cache(directory, signature):
            self.note('再利用 HISAT2インデックス')
            return prefix
        self.update(stage='hisat2_index', message='HISAT2のゲノム索引を作成中（初回のみ長時間かかります）', transfer=None, activity='')
        self.command(['hisat2-build', '-p', self.threads, fasta, prefix])
        paths = [Path(str(prefix) + f'.{i}.ht2') for i in range(1, 9)]
        if not all(p.is_file() for p in paths):
            paths = [Path(str(prefix) + f'.{i}.ht2l') for i in range(1, 9)]
        if not all(p.is_file() and p.stat().st_size for p in paths):
            raise RuntimeError('HISAT2インデックスの生成が不完全です。')
        mark_complete(directory, signature, paths)
        return prefix

    def check_bam(self, bam):
        self.command(['samtools', 'quickcheck', bam])
        if self.reference_fasta is not None:
            with rnako.pysam.FastaFile(str(self.reference_fasta)) as reference, rnako.pysam.AlignmentFile(str(bam), 'rb') as alignment:
                if dict(zip(reference.references, reference.lengths)) != dict(zip(alignment.references, alignment.lengths)):
                    raise RuntimeError(f'{bam.name}: BAMの参照名・全contig長が選択したFASTAと一致しません。再マッピングしてください。')

    def cleanup_fastqs(self, run):
        """Discard generated reads only after a usable run BAM is available."""
        from storage import Storage, is_link
        if not re.fullmatch(r'[SED]RR\d+', run):
            raise ValueError('Invalid run accession for FASTQ cleanup')
        directories = []
        raw = self.data / 'fastq' / run
        if raw.exists():
            Storage.safe_path(self.data, raw)
            directories.append(raw)
        trimmed = self.data / 'trimmed' / run
        if trimmed.exists():
            Storage.safe_path(self.data, trimmed)
            for directory in trimmed.iterdir():
                if directory.is_dir():
                    Storage.safe_path(self.data, directory)
                    directories.append(directory)
        targets = []
        markers = []
        for directory in directories:
            reads = [p for p in directory.iterdir() if p.name.startswith(run)
                     and p.name.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz'))]
            for path in reads:
                if is_link(path.lstat()) or not path.is_file():
                    raise ValueError('Unsafe FASTQ cleanup target')
            targets.extend(reads)
            marker = directory / 'complete.json'
            if reads and marker.exists():
                if is_link(marker.lstat()) or not marker.is_file():
                    raise ValueError('Unsafe FASTQ completion marker')
                markers.append(marker)
        # Invalidate caches first, including when deletion is interrupted.
        for marker in markers:
            marker.unlink()
        size = sum(p.stat().st_size for p in targets)
        for path in targets:
            path.unlink()
        if targets:
            self.note(f'{run}: BAM検証済みのため中間FASTQ {len(targets)}件を削除 ({size / 1024**3:.2f} GiB)。SRA・BAM・fastpレポートは保持。')
        return size

    def mapped_run(self, row, index, versions):
        run = row['run_accession']
        with self.exclusive('run:' + run, f'{run}: 同時に実行中の別の解析がこのrunを処理中のため、完了を待って再利用します。'):
            return self.locked_mapped_run(row, index, versions)

    def locked_mapped_run(self, row, index, versions):
        run = row['run_accession']
        self.current_run = run
        signature = {'version': PIPELINE_VERSION, 'run': run, 'layout': row['library_layout'],
                     'index': read_json(index.parent / 'complete.json'), 'hisat2': versions['hisat2'],
                     'samtools': versions['samtools'], 'options': mapping_options(row['library_layout']),
                     'chunk_reads': MAPPING_CHUNK_READS,
                     'preprocessing': {'fastp': versions['fastp'], 'fastp_version': FASTP_VERSION,
                                       'fastp_options': fastp_options(row['library_layout']),
                                       'fasterq': versions['fasterq-dump'], 'fastq_version': FASTQ_VERSION}}
        directory = self.data / 'bam' / run / digest(signature)
        directory.mkdir(parents=True, exist_ok=True)
        bam = directory / (run + '.bam')
        if valid_cache(directory, signature):
            self.check_bam(bam)
            self.preprocessing[run] = read_json(directory / 'preprocessing.json') or {}
            self.note(f'再利用 BAM: {run}')
            self.cleanup_fastqs(run)
            return bam
        if self.ram_root is None:
            return self.map_reads(row, index, directory, bam, signature, versions)
        # Download, convert, trim and map in RAM; only the finished BAM is written to disk.
        self.work_dir = self.ram_root / run
        shutil.rmtree(self.work_dir, ignore_errors=True)
        self.work_dir.mkdir(parents=True)
        try:
            return self.map_reads(row, index, directory, bam, signature, versions)
        finally:
            shutil.rmtree(self.work_dir, ignore_errors=True)
            self.work_dir = None
            self.release_ram()

    def map_reads(self, row, index, directory, bam, signature, versions):
        run = row['run_accession']
        self.update(stage='download', message=f'{run}: SRA取得・FASTQ変換' + ('（RAM上）' if self.work_dir else ''))
        raw = self.fastqs(row, versions)
        fastqs = self.trimmed_fastqs(row, raw, versions)
        if self.work_dir:
            for path in raw:
                path.unlink()  # free RAM before mapping; fastp output is all HISAT2 needs
            self.shrink_ram(sum(p.stat().st_size for p in fastqs) * RAM_PER_TRIMMED_BYTE + RAM_MAPPING_MARGIN)
        self.update(stage='mapping', message=f'{run}: HISAT2マッピング・BAM作成' + ('（RAM上）' if self.work_dir else ''))
        scratch = self.work / 'bam' / run if self.work_dir else directory
        scratch.mkdir(parents=True, exist_ok=True)
        partial = scratch / 'mapping.partial.bam'
        for stale in [*scratch.glob('sort.partial*'), *scratch.glob('chunk*.fq')]:
            stale.unlink()
        memory = available_memory()
        workers = self.shared.workers
        sort_threads, sort_memory = sort_resources(self.threads, memory and memory // workers)
        sort = ['samtools', 'sort', '-@', str(sort_threads), '-m', f'{sort_memory}M',
                '-T', str(scratch / 'sort.partial'), '-o', str(partial), '-']
        self.checkpoint()
        chunks = ReadChunks(fastqs, scratch)
        limit = hisat2_memory_limit(memory, workers)
        align = (['prlimit', f'--as={limit}', '--'] if limit else []) + [
            'hisat2', '-p', str(self.threads), *mapping_options(row['library_layout']),
            '-x', str(index), '--rg-id', run, '--rg', 'SM:' + run]
        inputs = lambda paths: ['-1', str(paths[0]), '-2', str(paths[1])] if len(paths) == 2 else ['-U', str(paths[0])]
        self.note('$ ' + ' '.join(align + inputs(chunks.sets[0])) + ' | ' + ' '.join(sort))
        self.note(f'{run}: {MAPPING_CHUNK_READS:,}リードごとにHISAT2を起動し直してマッピング（-k {MULTIMAP_MAX}で'
                  'HISAT2のメモリが増え続けるため）。入力: ' + ', '.join(map(str, fastqs)))
        with self.log.open('ab') as log:
            sorter = subprocess.Popen(sort, stdin=subprocess.PIPE, stdout=log, stderr=log,
                                      start_new_session=(os.name == 'posix'))
            aligner, block = None, 0
            try:
                while (paths := chunks.next()) is not None:
                    block += 1
                    # Only the first block writes the SAM header; samtools sort reads one continuous stream.
                    aligner = subprocess.Popen(align + inputs(paths) + ([] if block == 1 else ['--no-hd']),
                                               stdout=sorter.stdin, stderr=log, stdin=subprocess.DEVNULL,
                                               start_new_session=(os.name == 'posix'))
                    align_code = self.wait_process(aligner, 'HISAT2 / samtools sort')
                    aligner = None
                    if align_code:
                        if sorter.poll() is not None:
                            hint = f'samtools sortが先に終了しました（終了コード {sorter.returncode}）。'
                        elif limit:
                            hint = f'HISAT2のメモリ上限 {limit / GIB:.0f} GiB を超えた可能性があります。同時実行数を減らしてください。'
                        else:
                            hint = ''
                        raise RuntimeError(f'{run}: マッピング失敗（HISAT2 終了コード {align_code}、'
                                           f'{block}番目のブロック）。{hint}')
                    self.note(f'{run}: HISAT2 {block}番目のブロック完了（累計 {chunks.reads:,}リード）')
                sorter.stdin.close()
                sort_code = self.wait_process(sorter, 'HISAT2 / samtools sort')
            except BaseException:
                if aligner is not None:
                    stop_process(aligner)
                stop_process(sorter)
                raise
            finally:
                chunks.close()
                if not sorter.stdin.closed:
                    sorter.stdin.close()
        if sort_code:
            raise RuntimeError(f'{run}: マッピング失敗 (samtools sort={sort_code})')
        self.check_bam(partial)
        self.command(['samtools', 'index', '-@', self.threads, partial])
        if self.work_dir:
            # Keep the fastp QC reports beside the BAM; the RAM copies are discarded with the run.
            report = dict(self.preprocessing[run])
            for key in ('fastp_html', 'fastp_json'):
                kept = directory / Path(report[key]).name
                shutil.copyfile(report[key], kept)
                report[key] = str(kept)
            report.pop('raw_fastq', None)
            report.pop('trimmed_fastq', None)
            self.preprocessing[run] = report
            # Copy to disk under the partial name first so an interruption never leaves a truncated BAM.
            self.update(message=f'{run}: BAMをディスクへ保存中')
            for source in (partial, Path(str(partial) + '.bai')):
                shutil.copyfile(source, directory / source.name)
                source.unlink()
            partial = directory / partial.name
        os.replace(partial, bam)
        os.replace(Path(str(partial) + '.bai'), Path(str(bam) + '.bai'))
        preprocessing = directory / 'preprocessing.json'
        atomic_json(preprocessing, self.preprocessing[run])
        mark_complete(directory, signature, [bam, Path(str(bam) + '.bai'), preprocessing])
        self.command(['samtools', 'flagstat', bam])
        self.cleanup_fastqs(run)
        return bam

    def sample_bam(self, accession, bams, versions):
        if len(bams) == 1:
            return bams[0]
        with self.exclusive('sample:' + accession, f'{accession}: 同時に実行中の別の解析がこの検体を統合中のため待機しています。'):
            return self.locked_sample_bam(accession, bams, versions)

    def locked_sample_bam(self, accession, bams, versions):
        signature = {'version': PIPELINE_VERSION, 'sample': accession, 'runs': fingerprint(bams),
                     'samtools': versions['samtools'], 'merge': 'prefix-query-name-with-run'}
        directory = self.data / 'samples' / accession / digest(signature)
        directory.mkdir(parents=True, exist_ok=True)
        bam = directory / (accession + '.bam')
        if not valid_cache(directory, signature):
            # Read names can repeat between runs. Prefix before merging so gene/junction
            # fragment counting in rnako never collapses different runs' fragments.
            renamed = []
            for i, source in enumerate(bams):
                target = directory / f'run_{i}.bam'
                with rnako.pysam.AlignmentFile(str(source), 'rb') as src:
                    with rnako.pysam.AlignmentFile(str(target), 'wb', template=src) as dst:
                        for number, read in enumerate(src):
                            if number % 10000 == 0:
                                self.checkpoint()
                            name = f'{i}:{read.query_name}'
                            if len(name) > 254:
                                name = f'{i}:' + hashlib.sha256(read.query_name.encode()).hexdigest()
                            read.query_name = name
                            dst.write(read)
                renamed.append(target)
            partial = directory / 'merged.partial.bam'
            self.command(['samtools', 'merge', '-f', '-@', self.threads, partial, *renamed])
            self.command(['samtools', 'quickcheck', partial])
            self.command(['samtools', 'index', '-@', self.threads, partial])
            os.replace(partial, bam)
            os.replace(Path(str(partial) + '.bai'), Path(str(bam) + '.bai'))
            mark_complete(directory, signature, [bam, Path(str(bam) + '.bai')])
            for path in renamed:
                path.unlink()
        else:
            self.command(['samtools', 'quickcheck', bam])
            self.note(f'再利用 SRX統合BAM: {accession}')
        return bam

    def run(self, request, job_dir):
        self.update(stage='preflight', message='解析環境とSRX情報を確認中', progress=2)
        self.note(f'実行設定: CPUスレッド数 {self.threads}' + (' / RAM上で取得〜マッピング' if request.get('ram_work') else ''))
        if request.get('ram_work'):
            if self.ram is None or not os.path.ismount(self.ram):
                raise RuntimeError('RAM作業領域（tmpfs）がありません。start.shで起動し直すか、RAMを使わない設定で実行してください。')
            self.ram_root = self.ram / job_dir.name
            shutil.rmtree(self.ram_root, ignore_errors=True)  # left over from an interrupted attempt
            self.note(f'RAM作業領域: {self.ram}（{shutil.disk_usage(self.ram).total / GIB:.0f} GiB）。'
                      'SRA・FASTQはRAM上で処理し、BAMとfastpレポートだけをディスクに保存します。')
        versions = self.tools()
        rows = []
        seen_runs = set()
        for label in ('wt', 'ko'):
            for accession in request[label]:
                self.update(message=f'{accession}: SRXとrunの対応を確認中')
                runs = resolve_runs(accession, request['assembly'], self.data, self.note)
                if any(r['run_accession'] in seen_runs for r in runs):
                    raise ValueError('複数のSRXに同じrunが含まれています。')
                seen_runs.update(r['run_accession'] for r in runs)
                rows.append((accession, label.upper(), runs))
                self.note(f'{accession} ({label.upper()}): ' + ', '.join(r['run_accession'] for r in runs))
        self.update(stage='reference', message='参照FASTA・GTFの取得／再利用、HISAT2索引を準備中', progress=5)
        waiting = '同時に実行中の別の解析が参照・索引を準備中のため、完了を待って再利用します。'
        with self.exclusive('reference', waiting):
            fasta, gtf = rnako.prepare_reference(request['assembly'], self.references,
                                                progress=self.update, check_pause=self.checkpoint)
        self.update(stage='reference_annotation', message='対象遺伝子の注釈を確認中', transfer=None)
        rnako.load_gene(gtf, request['gene'])  # Fail before downloading runs for a misspelled gene.
        with self.exclusive('reference', waiting):
            index = self.index(fasta, versions)
        self.update(reference_ready=True, activity='参照ゲノムとHISAT2索引を準備できました。')
        if request.get('pause_after_reference'):
            raise PipelinePaused('参照・索引の準備が完了したため中断しました。「再開」でデータ取得へ進みます。')
        samples = []
        done = 0
        try:
            for accession, label, runs in rows:
                bams = []
                for row in runs:
                    self.update(progress=10 + int(65 * done / len(seen_runs)), runs_done=done,
                                runs_total=len(seen_runs), message=f'{accession} / {row["run_accession"]} を処理中')
                    bams.append(self.mapped_run(row, index, versions))
                    done += 1
                    self.update(runs_done=done)
                self.update(stage='merge', message=f'{accession}: runを検体単位に統合中')
                samples.append({'srx_id': accession, 'metadata_label': label,
                                'bam_path': str(self.sample_bam(accession, bams, versions))})
        finally:
            if self.ram_root is not None:
                shutil.rmtree(self.ram_root, ignore_errors=True)
        input_path = job_dir / 'input_samples.csv'
        with input_path.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['srx_id', 'metadata_label', 'bam_path'])
            writer.writeheader()
            writer.writerows(samples)
        atomic_json(job_dir / 'provenance.json', {'request': request, 'tools': versions, 'pipeline_version': PIPELINE_VERSION,
                    'threads': self.threads, 'ram_work': self.ram_root is not None,
                    'experiments': [{'accession': a, 'label': l, 'runs': r} for a, l, r in rows],
                    'reference': fingerprint([fasta, gtf]), 'hisat2_index': str(index),
                    'workflow': ['prefetch', 'vdb-validate', 'fasterq-dump --split-files', 'fastp',
                                 'hisat2', 'samtools sort', 'samtools quickcheck', 'samtools index'],
                    'alignment_options': '-k 100 --secondary; paired: --no-mixed --no-discordant; unstranded',
                    'preprocessing': self.preprocessing,
                    'aggregation': 'all public runs per experiment merged; run-prefixed read names'})
        self.update(stage='analysis', message='既存のKO解析を実行中（全BAMのライブラリ計数を含みます）', progress=80)
        self.command([sys.executable, ROOT / 'rnako.py', '--gene', request['gene'], '--assembly', request['assembly'],
                      '--samples', input_path, '--fasta', fasta, '--gtf', gtf, '--out', job_dir / 'output', '--skip-region-analysis']
                     + (['--confirm-replicates'] if request.get('confirm_replicates') is True else []))
        self.update(stage='analysis', message='候補変異を全検体で再観測し、比較表を作成中', progress=95)
        write_observations(samples, job_dir / 'output', fasta, gtf, request['gene'], self.checkpoint)
        self.update(message='エクソンごとのカバレッジとWT比較を作成中', progress=98)
        write_exon_coverage(samples, job_dir / 'output', fasta, gtf, request['gene'], self.checkpoint)
        self.update(message='RNA領域欠損と条件付きアミノ酸配列を評価中', progress=98)
        from region_loss import write_region_loss
        write_region_loss(job_dir / 'output', fasta, gtf, request['gene'], self.checkpoint)
        self.checkpoint()
        write_conclusion(job_dir / 'output', request.get('confirm_replicates') is True)
        self.update(stage='completed', message='解析が完了しました', progress=100, status='completed')


def write_observations(samples, output, fasta, gtf, gene, check_pause=None):
    """Expose the core's candidate-site recheck for every sample, including 0 alt."""
    with (output / 'events.csv').open(newline='', encoding='utf-8') as handle:
        events = {r['event_id']: r for r in csv.DictReader(handle)}
    columns = ['srx_id', 'metadata_label', 'event_id', 'chrom', 'pos_1based', 'ref', 'alt',
               'depth_at_anchor', 'alt_reads', 'local_alt_fraction', 'coverage_status']
    transcripts = rnako.load_gene(gtf, gene)
    annotated_chrom = next(iter(transcripts.values())).chrom
    settings = read_json(output / 'run.json')['thresholds']
    with rnako.pysam.FastaFile(str(fasta)) as fa, (output / 'observations.csv').open('w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        for sample in samples:
            if check_pause:
                check_pause()
            with rnako.pysam.AlignmentFile(sample['bam_path'], 'rb') as bam:
                chrom = rnako.resolve_contig(annotated_chrom, bam, fa)
                for eid, event in events.items():
                    depth, alt = rnako.observe_event(bam, chrom, event, settings['min_mapq'], settings['min_bq'])
                    writer.writerow({'srx_id': sample['srx_id'], 'metadata_label': sample['metadata_label'],
                        'event_id': eid, **{k: event[k] for k in ('chrom', 'pos_1based', 'ref', 'alt')},
                        'depth_at_anchor': depth, 'alt_reads': alt,
                        'local_alt_fraction': round(alt / depth, 5) if depth else '',
                        'coverage_status': 'covered' if depth >= settings['min_event_depth'] else 'low_depth'})


def write_exon_coverage(samples, output, fasta, gtf, gene, check_pause=None):
    """Descriptive RNA coverage evidence; never infer a genomic exon deletion."""
    check_pause = check_pause or (lambda: None)
    check_pause()
    from collections import Counter, defaultdict
    from statistics import mean

    transcripts = rnako.load_gene(gtf, gene)
    annotated_chrom = next(iter(transcripts.values())).chrom
    annotations = defaultdict(list)
    for transcript in transcripts.values():
        for number, interval in enumerate(transcript.ordered_exons(), 1):
            annotations[interval].append(f'{transcript.tid}:exon{number}')
    intervals = sorted(annotations)
    merged = rnako.gene_exons(transcripts)
    gene_length = sum(e - s for s, e in merged)
    settings = read_json(output / 'run.json')['thresholds']
    min_depth = settings['min_event_depth']
    rows = []
    with rnako.pysam.FastaFile(str(fasta)) as fa:
        for sample in samples:
            check_pause()
            depths = Counter()
            junctions = defaultdict(set)
            fragments = defaultdict(set)
            with rnako.pysam.AlignmentFile(sample['bam_path'], 'rb') as bam:
                chrom = rnako.resolve_contig(annotated_chrom, bam, fa)
                for col in bam.pileup(chrom, merged[0][0], merged[-1][1], truncate=True,
                                      stepper='nofilter', min_base_quality=0, max_depth=1000000,
                                      ignore_overlaps=True):
                    check_pause()
                    pos = col.reference_pos
                    if not any(s <= pos < e for s, e in merged):
                        continue
                    for p in col.pileups:
                        r, q = p.alignment, p.query_position
                        if (not rnako.read_ok(r, settings['min_mapq']) or p.is_refskip or p.is_del
                                or q is None or r.query_qualities is None
                                or r.query_qualities[q] < settings['min_bq']
                                or r.query_sequence[q].upper() not in 'ACGT'):
                            continue
                        depths[pos] += 1
                for read in bam.fetch(chrom, merged[0][0], merged[-1][1]):
                    if not rnako.read_ok(read, settings['min_mapq']):
                        continue
                    for interval in intervals:
                        if rnako.exonic_overlap(read, [interval]):
                            fragments[interval].add(read.query_name)
                    for junction in rnako.junctions_from_read(read):
                        junctions[junction].add(read.query_name)
            total_depth = sum(depths.values())
            total_sufficient = sum(d >= min_depth for d in depths.values())
            for s, e in intervals:
                length = e - s
                values = [depths[pos] for pos in range(s, e)]
                depth_sum = sum(values)
                sufficient = sum(d >= min_depth for d in values)
                rest_length = gene_length - length
                rest_mean = (total_depth - depth_sum) / rest_length if rest_length else None
                rest_fraction = (total_sufficient - sufficient) / rest_length if rest_length else 0
                exon_mean = depth_sum / length
                relative = exon_mean / rest_mean if rest_mean else None
                longest_start = longest_length = current = 0
                for pos, d in enumerate(values):
                    current = 0 if d else current + 1
                    if current > longest_length:
                        longest_length = current
                        longest_start = s + pos - current + 2
                skip_names = set()
                for (a, b), names in junctions.items():
                    if a <= s and b >= e:
                        skip_names.update(names)
                rows.append({
                    'srx_id': sample['srx_id'], 'metadata_label': sample['metadata_label'],
                    'gene': next(iter(transcripts.values())).gene_name, 'chrom': annotated_chrom,
                    'exon_start_1based': s + 1, 'exon_end_1based': e,
                    'transcript_exons': ';'.join(sorted(annotations[(s, e)])),
                    'exon_length': length, 'overlapping_fragments': len(fragments[(s, e)]),
                    'mean_depth': round(exon_mean, 3), 'min_depth': min(values),
                    'covered_fraction': round(sum(d > 0 for d in values) / length, 5),
                    'sufficient_depth_fraction': round(sufficient / length, 5),
                    'longest_zero_depth_region': (f'{annotated_chrom}:{longest_start}-'
                                                  f'{longest_start + longest_length - 1}') if longest_length else '',
                    'skipping_junction_fragments': len(skip_names),
                    'other_exons_mean_depth': round(rest_mean, 3) if rest_mean is not None else '',
                    'other_exons_sufficient_fraction': round(rest_fraction, 5),
                    'relative_exon_depth': relative,
                })
    for row in rows:
        reference = [r['relative_exon_depth'] for r in rows
                     if r['metadata_label'] == 'WT' and r['srx_id'] != row['srx_id']
                     and r['exon_start_1based'] == row['exon_start_1based']
                     and r['exon_end_1based'] == row['exon_end_1based']
                     and r['sufficient_depth_fraction'] >= .8
                     and r['other_exons_sufficient_fraction'] >= .8
                     and r['relative_exon_depth'] is not None]
        baseline = mean(reference) if reference else None
        ratio = row['relative_exon_depth'] / baseline if baseline and row['relative_exon_depth'] is not None else None
        row['wt_reference_samples'] = len(reference)
        row['relative_depth_vs_wt'] = round(ratio, 5) if ratio is not None else ''
        row['coverage_status'] = ('no_coverage' if row['covered_fraction'] == 0 else
                                  'low_coverage' if row['sufficient_depth_fraction'] < .8 else 'covered')
        if row['metadata_label'] == 'WT':
            status = 'reference_sample'
        elif baseline is None:
            status = 'insufficient_wt_coverage'
        elif row['other_exons_sufficient_fraction'] < .8 or ratio is None:
            status = 'insufficient_sample_coverage'
        elif ratio <= .2 and row['covered_fraction'] < .5:
            status = 'suspected_local_depletion'
        elif row['sufficient_depth_fraction'] < .8:
            status = 'insufficient_sample_coverage'
        else:
            status = 'no_specific_depletion'
        row['comparison_status'] = status
    for row in rows:
        row['relative_exon_depth'] = round(row['relative_exon_depth'], 5) if row['relative_exon_depth'] is not None else ''
    target = output / 'exon_coverage.csv'
    temporary = target.with_suffix('.csv.partial')
    with temporary.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, target)

