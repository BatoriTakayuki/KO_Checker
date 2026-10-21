"""Manage only generated cache directories, with preview-bound deletion."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import re
import secrets
import shutil
import stat
import time


KINDS = {
    'sra': ('SRA', '再利用時はSRAの再取得が必要です。'),
    'fastq': ('元FASTQ', '再利用時はSRAからFASTQへ再変換します。'),
    'trimmed': ('前処理済みFASTQ・QC', '再利用時はfastpの前処理をやり直します。'),
    'bam': ('runごとのBAM', '再利用時はマッピングをやり直します。'),
    'samples': ('検体の統合BAM', '再利用時はrunのBAMから統合し直します。'),
    'reference': ('参照FASTA・注釈', '再利用時は参照の取得・展開が必要です。'),
    'index': ('HISAT2索引', '再利用時は索引を構築し直します。'),
}
RUN = re.compile(r'(?:SRR|ERR|DRR)[0-9]+\Z')
EXPERIMENT = re.compile(r'(?:SRX|ERX|DRX)[0-9]+\Z')
INDEX = re.compile(r'[a-f0-9]{24}\Z')
MAX_SELECTION = 300


def is_link(info):
    # Junctions and other Windows reparse points must never be traversed.
    return stat.S_ISLNK(info.st_mode) or bool(getattr(info, 'st_file_attributes', 0) & 0x400)


class StorageBusy(ValueError):
    pass


class Storage:
    def __init__(self, data, references):
        self.data = Path(os.path.abspath(data))
        self.references = Path(os.path.abspath(references))
        self.previews = {}

    @staticmethod
    def safe_path(root, path):
        """Check absolute containment and each ancestor before recursive work."""
        if path == root or not path.is_relative_to(root):
            raise ValueError('保存先の外側は削除できません。')
        current = root
        for component in (None, *path.relative_to(root).parts):
            if component is not None:
                current = current / component
            if is_link(current.lstat()):
                raise ValueError('リンク・ジャンクションを含むデータは画面から削除できません。')
        resolved_root, resolved = root.resolve(strict=True), path.resolve(strict=True)
        if resolved == resolved_root or not resolved.is_relative_to(resolved_root) or os.path.ismount(path):
            raise ValueError('保存先の外側やマウント先は削除できません。')
        if not stat.S_ISDIR(path.lstat().st_mode):
            raise ValueError('削除対象はキャッシュのフォルダに限ります。')
        return resolved

    @staticmethod
    def scan(path):
        size = files = 0
        digest = hashlib.sha256()
        pending = [path]
        while pending:
            current = pending.pop()
            info = current.lstat()
            if is_link(info) or (current != path and os.path.ismount(current)):
                raise ValueError('リンク・ジャンクション・別のマウント先を含むため削除できません。')
            digest.update(repr((str(current.relative_to(path)), info.st_dev, info.st_ino,
                                info.st_mode, info.st_size, info.st_mtime_ns)).encode('utf-8'))
            if stat.S_ISDIR(info.st_mode):
                with os.scandir(current) as entries:
                    pending.extend(sorted((Path(e.path) for e in entries), reverse=True))
            elif stat.S_ISREG(info.st_mode):
                size += info.st_size
                files += 1
            else:
                raise ValueError('通常のファイル以外を含むため削除できません。')
        return size, files, digest.hexdigest()

    def candidates(self):
        for kind in ('sra', 'fastq', 'trimmed', 'bam', 'samples'):
            parent = self.data / kind
            if not parent.exists() and not parent.is_symlink():
                continue
            self.safe_path(self.data, parent)
            pattern = EXPERIMENT if kind == 'samples' else RUN
            for path in sorted(parent.iterdir()):
                if pattern.fullmatch(path.name):
                    yield kind, self.data, path, path.name
        for assembly in ('mm10', 'mm39', 'hg38'):
            path = self.references / assembly
            if path.exists() or path.is_symlink():
                yield 'reference', self.references, path, assembly
        parent = self.references / 'hisat2'
        if parent.exists() or parent.is_symlink():
            self.safe_path(self.references, parent)
            for path in sorted(parent.iterdir()):
                if INDEX.fullmatch(path.name):
                    yield 'index', self.references, path, path.name

    def catalog(self):
        entries, warnings = [], []
        try:
            for kind, root, path, label in self.candidates():
                relative = str(path.relative_to(root)).replace('\\', '/')
                entry = {'id': hashlib.sha256(f'{kind}:{relative}'.encode()).hexdigest()[:24],
                         'kind': kind, 'kind_label': KINDS[kind][0], 'label': label,
                         'path': ('references/' if kind in ('reference', 'index') else 'workdata/') + relative,
                         'description': KINDS[kind][1], 'deletable': True,
                         'size_bytes': 0, 'file_count': 0}
                try:
                    self.safe_path(root, path)
                    entry['size_bytes'], entry['file_count'], signature = self.scan(path)
                    entry['_signature'] = signature
                except (OSError, ValueError) as error:
                    entry.update(deletable=False, error=str(error))
                entry['_root'], entry['_path'] = root, path
                entries.append(entry)
        except (OSError, ValueError) as error:
            warnings.append(str(error))
        return entries, warnings

    @staticmethod
    def public(entry):
        return {k: v for k, v in entry.items() if not k.startswith('_')}

    def inventory(self):
        entries, warnings = self.catalog()
        totals = {kind: sum(e['size_bytes'] for e in entries if e['kind'] == kind) for kind in KINDS}
        disk = self.data
        while not disk.exists() and disk != disk.parent:
            disk = disk.parent
        return {'entries': [self.public(e) for e in entries], 'totals': totals,
                'total_bytes': sum(totals.values()), 'free_bytes': shutil.disk_usage(disk).free,
                'warnings': warnings, 'max_selection': MAX_SELECTION}

    def preview(self, ids):
        if (not isinstance(ids, list) or not 1 <= len(ids) <= MAX_SELECTION
                or any(not isinstance(i, str) or not re.fullmatch(r'[a-f0-9]{24}', i) for i in ids)
                or len(ids) != len(set(ids))):
            raise ValueError(f'削除対象を重複なく1〜{MAX_SELECTION}件選択してください。')
        entries, _ = self.catalog()
        available = {e['id']: e for e in entries}
        selected = []
        for key in ids:
            entry = available.get(key)
            if entry is None or not entry['deletable']:
                raise ValueError('削除対象が変更されたか、安全に削除できません。容量一覧を更新してください。')
            selected.append(entry)
        current = time.monotonic()
        self.previews = {k: v for k, v in self.previews.items() if v['expires'] > current}
        if len(self.previews) >= 16:
            self.previews.pop(next(iter(self.previews)))
        token = secrets.token_urlsafe(32)
        self.previews[token] = {'expires': current + 300, 'entries': selected}
        return {'token': token, 'entries': [self.public(e) for e in selected],
                'total_bytes': sum(e['size_bytes'] for e in selected), 'expires_in_seconds': 300}

    def delete(self, token, confirmed):
        if confirmed is not True or not isinstance(token, str):
            raise ValueError('削除内容を確認してから実行してください。')
        preview = self.previews.pop(token, None)
        if preview is None or preview['expires'] <= time.monotonic():
            raise ValueError('削除の確認が期限切れです。対象と容量を確認し直してください。')
        # Validate every target before deleting any of them. Jobs holds its lock
        # throughout this call, preventing submissions/resumes from racing cleanup.
        selected = preview['entries']
        for entry in selected:
            try:
                self.safe_path(entry['_root'], entry['_path'])
                _, _, signature = self.scan(entry['_path'])
                if signature != entry['_signature']:
                    raise ValueError('対象の内容が変更されました。容量一覧を更新して確認し直してください。')
            except OSError:
                raise ValueError('対象が変更されました。容量一覧を更新して確認し直してください。') from None
        deleted, failed = [], []
        for entry in selected:
            try:
                # Recheck the resolved absolute target immediately before rmtree.
                self.safe_path(entry['_root'], entry['_path'])
                shutil.rmtree(entry['_path'])
                deleted.append(self.public(entry))
            except (OSError, ValueError) as error:
                failed.append({**self.public(entry), 'error': str(error)})
        return {'deleted': deleted, 'failed': failed,
                'deleted_bytes': sum(e['size_bytes'] for e in deleted)}
