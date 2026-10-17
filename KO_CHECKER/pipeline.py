"""SRX -> SRA runs -> FASTQ -> HISAT2 BAM -> existing RNA evidence analysis."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import rnako

PIPELINE_VERSION = 1
ACCESSION = re.compile(r"(?:SRX|ERX|DRX)[0-9]+\Z")
RUN = re.compile(r"(?:SRR|ERR|DRR)[0-9]+\Z")
ROOT = Path(__file__).resolve().parent


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
    temporary = path.with_suffix(path.suffix + '.partial')
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
    def __init__(self, data, references, threads, update, log, should_pause=None):
        self.data = Path(data).resolve()
        self.references = Path(references).resolve()
        self.threads = threads
        self._update = update
        self.log = Path(log)
        self.should_pause = should_pause or (lambda: False)
        self.last_heartbeat = 0
        self.current_run = None
        self.current_index = None

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
            folder = self.data / 'sra' / self.current_run
            paths = (p for p in folder.rglob('*') if p.is_file() and p.suffix not in ('.lock', '.prf'))
            storage_label = 'SRA保存領域'
        elif command == 'fasterq-dump' and self.current_run:
            paths = (self.data / 'fastq' / self.current_run).glob('*.fastq')
            storage_label = 'FASTQ書き出し済み'
        elif command == 'HISAT2 / samtools sort' and self.current_run:
            paths = (self.data / 'bam' / self.current_run).rglob('*.partial.bam')
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
        names = ('prefetch', 'vdb-validate', 'fasterq-dump', 'hisat2', 'hisat2-build', 'samtools')
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
        signature = {'version': PIPELINE_VERSION, 'run': run, 'layout': row['library_layout'],
                     'fasterq': versions['fasterq-dump'], 'mode': 'split-files'}
        directory = self.data / 'fastq' / run
        directory.mkdir(parents=True, exist_ok=True)
        paired = row['library_layout'] == 'PAIRED'
        paths = [directory / (run + suffix + '.fastq') for suffix in (('_1', '_2') if paired else ('',))]
        if valid_cache(directory, signature):
            self.note(f'再利用 FASTQ: {run}')
            return paths
        archive_root = self.data / 'sra'
        archive_root.mkdir(parents=True, exist_ok=True)
        # prefetch resumes interrupted downloads and validates archive availability.
        self.update(stage='sra_download', message=f'{run}: SRAデータを取得中', transfer=None)
        self.command(['prefetch', run, '--max-size', 'u', '-O', archive_root / run])
        # Depending on prefetch's output-directory handling, -O <root>/<run>
        # places the archive directly there or creates another <run> directory.
        # Only accept this accession's exact filename in these two known layouts.
        candidates = [parent / (run + suffix)
                      for parent in (archive_root / run, archive_root / run / run)
                      for suffix in ('.sra', '.sralite')]
        archive = next((p for p in candidates if p.is_file() and p.stat().st_size), None)
        if archive is None:
            raise RuntimeError(f'{run}: prefetch後のSRAファイルが見つかりません。確認先: '
                               + ', '.join(map(str, candidates)))
        self.update(stage='sra_validate', message=f'{run}: 取得データの整合性を検査中')
        self.command(['vdb-validate', archive])
        temp = directory / 'scratch'
        temp.mkdir(exist_ok=True)
        self.update(stage='fastq_convert', message=f'{run}: SRAをFASTQへ変換中')
        self.command(['fasterq-dump', archive, '--split-files', '--force', '-e', self.threads,
                      '-O', directory, '-t', temp])
        if not paired and not paths[0].exists() and (directory / (run + '_1.fastq')).exists():
            os.replace(directory / (run + '_1.fastq'), paths[0])
        if not all(p.is_file() and p.stat().st_size for p in paths):
            raise RuntimeError(f'{run}: 期待するFASTQが揃いません。ライブラリ構成を確認してください。')
        unexpected = [p for p in directory.glob('*.fastq') if p not in paths and p.stat().st_size]
        if unexpected:
            raise RuntimeError(f'{run}: 追加のFASTQが存在します。特殊なライブラリ構成には未対応です。')
        mark_complete(directory, signature, paths)
        return paths

    def index(self, fasta, versions):
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

    def mapped_run(self, row, index, versions):
        run = row['run_accession']
        self.current_run = run
        signature = {'version': PIPELINE_VERSION, 'run': run, 'layout': row['library_layout'],
                     'index': read_json(index.parent / 'complete.json'), 'hisat2': versions['hisat2'],
                     'samtools': versions['samtools'], 'options': ['--dta', '--rg-id', run]}
        directory = self.data / 'bam' / run / digest(signature)
        directory.mkdir(parents=True, exist_ok=True)
        bam = directory / (run + '.bam')
        if valid_cache(directory, signature):
            self.command(['samtools', 'quickcheck', bam])
            self.note(f'再利用 BAM: {run}')
            return bam
        self.update(stage='download', message=f'{run}: SRA取得・FASTQ変換')
        fastqs = self.fastqs(row, versions)
        self.update(stage='mapping', message=f'{run}: HISAT2マッピング・BAM作成')
        partial = directory / 'mapping.partial.bam'
        align = ['hisat2', '-p', str(self.threads), '--dta', '-x', str(index), '--rg-id', run, '--rg', 'SM:' + run]
        align += ['-1', str(fastqs[0]), '-2', str(fastqs[1])] if len(fastqs) == 2 else ['-U', str(fastqs[0])]
        sort = ['samtools', 'sort', '-@', str(self.threads), '-o', str(partial), '-']
        self.note('$ ' + ' '.join(align) + ' | ' + ' '.join(sort))
        self.checkpoint()
        with self.log.open('ab') as log:
            aligner = subprocess.Popen(align, stdout=subprocess.PIPE, stderr=log, stdin=subprocess.DEVNULL,
                                       start_new_session=(os.name == 'posix'))
            sorter = None
            try:
                sorter = subprocess.Popen(sort, stdin=aligner.stdout, stdout=log, stderr=log,
                                          start_new_session=(os.name == 'posix'))
                aligner.stdout.close()
                sort_code = self.wait_process(sorter, 'HISAT2 / samtools sort')
                align_code = self.wait_process(aligner, 'HISAT2')
            except BaseException:
                stop_process(aligner)
                if sorter is not None:
                    stop_process(sorter)
                raise
        if sort_code or align_code:
            raise RuntimeError(f'{run}: マッピング失敗 (HISAT2={align_code}, samtools={sort_code})')
        self.command(['samtools', 'quickcheck', partial])
        self.command(['samtools', 'index', '-@', self.threads, partial])
        os.replace(partial, bam)
        os.replace(Path(str(partial) + '.bai'), Path(str(bam) + '.bai'))
        mark_complete(directory, signature, [bam, Path(str(bam) + '.bai')])
        self.command(['samtools', 'flagstat', bam])
        return bam

    def sample_bam(self, accession, bams, versions):
        if len(bams) == 1:
            return bams[0]
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
        self.note(f'実行設定: CPUスレッド数 {self.threads}')
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
        fasta, gtf = rnako.prepare_reference(request['assembly'], self.references,
                                            progress=self.update, check_pause=self.checkpoint)
        self.update(stage='reference_annotation', message='対象遺伝子の注釈を確認中', transfer=None)
        rnako.load_gene(gtf, request['gene'])  # Fail before downloading runs for a misspelled gene.
        index = self.index(fasta, versions)
        self.update(reference_ready=True, activity='参照ゲノムとHISAT2索引を準備できました。')
        if request.get('pause_after_reference'):
            raise PipelinePaused('参照・索引の準備が完了したため中断しました。「再開」でデータ取得へ進みます。')
        samples = []
        done = 0
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
        input_path = job_dir / 'input_samples.csv'
        with input_path.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['srx_id', 'metadata_label', 'bam_path'])
            writer.writeheader()
            writer.writerows(samples)
        atomic_json(job_dir / 'provenance.json', {'request': request, 'tools': versions, 'pipeline_version': PIPELINE_VERSION,
                    'threads': self.threads,
                    'experiments': [{'accession': a, 'label': l, 'runs': r} for a, l, r in rows],
                    'reference': fingerprint([fasta, gtf]), 'hisat2_index': str(index),
                    'alignment_options': '--dta; unstranded; no adapter trimming',
                    'aggregation': 'all public runs per experiment merged; run-prefixed read names'})
        self.update(stage='analysis', message='既存のKO解析を実行中（全BAMのライブラリ計数を含みます）', progress=80)
        self.command([sys.executable, ROOT / 'rnako.py', '--gene', request['gene'], '--assembly', request['assembly'],
                      '--samples', input_path, '--fasta', fasta, '--gtf', gtf, '--out', job_dir / 'output'])
        self.update(stage='analysis', message='候補変異を全検体で再観測し、比較表を作成中', progress=95)
        write_observations(samples, job_dir / 'output', fasta, gtf, request['gene'], self.checkpoint)
        self.update(message='エクソンごとのカバレッジとWT比較を作成中', progress=98)
        write_exon_coverage(samples, job_dir / 'output', fasta, gtf, request['gene'], self.checkpoint)
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

