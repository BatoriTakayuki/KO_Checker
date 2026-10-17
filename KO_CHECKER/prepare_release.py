"""Build a source-only ZIP for GitHub using an explicit file allowlist."""
from __future__ import annotations

import hashlib
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent
PUBLISH_FILES = (
    'README.md', 'CHANGELOG.md', 'docs/analysis.md', 'docs/github-sharing.md',
    '.gitignore', '.dockerignore', '.gitattributes', '.env.example',
    '.github/workflows/tests.yml',
    'Dockerfile', 'compose.yaml', 'start.bat', 'start.ps1', 'requirements.txt', 'example_samples.csv',
    'rnako.py', 'pipeline.py', 'webapp.py', 'transfers.py',
    'static/index.html', 'static/style.css', 'static/app.js',
    'test_core.py', 'test_integration.py', 'test_web.py', 'test_resume.py', 'prepare_release.py',
)


def main():
    contents = {}
    for name in PUBLISH_FILES:
        path = ROOT / name
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise RuntimeError(f'Missing or unsafe source file: {name}')
        if path.stat().st_size > 2_000_000:
            raise RuntimeError(f'Unexpectedly large source file: {name}')
        contents[name] = path.read_bytes()
    # License selection belongs to the repository owner. Include one if added later.
    license_path = ROOT / 'LICENSE'
    if license_path.is_file() and not license_path.is_symlink():
        if license_path.stat().st_size > 100_000:
            raise RuntimeError('Unexpectedly large LICENSE')
        contents['LICENSE'] = license_path.read_bytes()
    manifest = ''.join(f'{hashlib.sha256(data).hexdigest()}  {name}\n'
                       for name, data in sorted(contents.items()))
    contents['SHA256SUMS.txt'] = manifest.encode('utf-8')
    manifest_path = ROOT / 'SHA256SUMS.txt'
    manifest_temporary = manifest_path.with_suffix('.txt.partial')
    manifest_temporary.write_bytes(contents['SHA256SUMS.txt'])
    manifest_temporary.replace(manifest_path)
    destination = ROOT / 'dist' / 'rnako-github.zip'
    destination.parent.mkdir(exist_ok=True)
    temporary = destination.with_suffix('.zip.partial')
    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(contents.items()):
            archive.writestr(name, data)
    temporary.replace(destination)
    print(f'Created: {destination}')
    print(f'Files: {len(contents)}, ZIP size: {destination.stat().st_size:,} bytes')


if __name__ == '__main__':
    main()
