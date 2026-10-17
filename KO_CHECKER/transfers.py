"""Resumable reference downloads with validators and byte-level progress."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import time
import urllib.error
import urllib.request


def download_file(url, destination, progress=None, check_pause=None):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.partial')
    metadata = partial.with_name(partial.name + '.json')
    check_pause = check_pause or (lambda: None)
    progress = progress or (lambda **kw: None)
    for attempt in range(3):
        check_pause()
        try:
            saved = json.loads(metadata.read_text(encoding='utf-8'))
        except (OSError, ValueError):
            saved = {}
        if not isinstance(saved, dict):
            saved = {}
        offset = partial.stat().st_size if partial.exists() else 0
        validator = saved.get('validator')
        if saved.get('url') != url or not validator:
            offset = 0
        headers = {'User-Agent': 'KO_Checker/0.2', 'Accept-Encoding': 'identity'}
        if offset:
            headers.update({'Range': f'bytes={offset}-', 'If-Range': validator})
        response = None
        try:
            response = urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=30)
            with response as source:
                status = source.status
                total = None
                if status == 206:
                    match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', source.headers.get('Content-Range', ''))
                    if (not offset or not match or int(match[1]) != offset
                            or int(match[2]) != int(match[3]) - 1
                            or (saved.get('total') is not None and int(match[3]) != saved['total'])):
                        # Never append a response whose range does not match our prefix.
                        metadata.unlink(missing_ok=True)
                        raise OSError('Invalid Content-Range; retrying from the beginning')
                    total = int(match[3])
                    returned_validator = source.headers.get(saved.get('validator_header', 'ETag'))
                    if returned_validator and returned_validator != validator:
                        metadata.unlink(missing_ok=True)
                        raise OSError('Reference changed; retrying from the beginning')
                elif status == 200:
                    offset = 0  # Range ignored or If-Range changed: replace, do not append.
                    length = source.headers.get('Content-Length')
                    total = int(length) if length else None
                else:
                    raise OSError(f'Unexpected HTTP status {status}')
                etag = source.headers.get('ETag', '')
                modified = source.headers.get('Last-Modified', '')
                if etag and not etag.startswith('W/'):
                    new_validator, validator_header = etag, 'ETag'
                elif modified:
                    new_validator, validator_header = modified, 'Last-Modified'
                elif offset:
                    new_validator, validator_header = validator, saved.get('validator_header', 'ETag')
                else:
                    new_validator, validator_header = None, None
                state = {'url': url, 'validator': new_validator,
                         'validator_header': validator_header, 'total': total}
                # Truncate before publishing a new validator, so a restart cannot combine
                # an old prefix with metadata describing a new remote representation.
                with partial.open('ab' if offset else 'wb') as target:
                    temporary = metadata.with_suffix('.json.partial')
                    temporary.write_text(json.dumps(state), encoding='utf-8')
                    os.replace(temporary, metadata)
                    done = offset
                    started = time.monotonic()
                    notified = 0.0
                    progress(name=destination.name, done_bytes=done, total_bytes=total,
                             rate_bps=0, resumed_from=offset)
                    while True:
                        check_pause()
                        read = getattr(source, 'read1', source.read)
                        chunk = read(256 * 1024)
                        if not chunk:
                            break
                        target.write(chunk)
                        target.flush()
                        done += len(chunk)
                        current = time.monotonic()
                        if current - notified >= 1:
                            progress(name=destination.name, done_bytes=done, total_bytes=total,
                                     rate_bps=(done - offset) / max(current - started, 0.001),
                                     resumed_from=offset)
                            notified = current
                    if not done or (total is not None and done != total):
                        raise OSError(f'Incomplete download: {done}/{total} bytes')
                    progress(name=destination.name, done_bytes=done, total_bytes=total,
                             rate_bps=(done - offset) / max(time.monotonic() - started, 0.001),
                             resumed_from=offset)
                check_pause()
                os.replace(partial, destination)
                metadata.unlink(missing_ok=True)
                return
        except urllib.error.HTTPError as error:
            if error.code == 416:
                metadata.unlink(missing_ok=True)  # Re-fetch the full representation safely.
            elif error.code not in (408, 429, 500, 502, 503, 504):
                raise
            if attempt == 2:
                raise
        except (OSError, urllib.error.URLError):
            if attempt == 2:
                raise
        # Partial data and its validator deliberately survive pause/network failure.
        for _ in range(10 * (attempt + 1)):
            check_pause()
            time.sleep(0.1)
