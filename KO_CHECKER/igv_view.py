"""Local IGV configuration; expose only files recorded by this analysis."""
import csv
import json
from pathlib import Path


def confined_file(path, roots):
    target = Path(path).resolve()
    if not any(target.is_relative_to(Path(root).resolve()) for root in roots):
        raise ValueError('File is outside analysis storage')
    if not target.is_file():
        raise FileNotFoundError('解析ファイルがありません。キャッシュ削除済みの場合は再解析してください。')
    return target


def assets(output, data, references):
    output = Path(output)
    settings = json.loads((output / 'run.json').read_text(encoding='utf-8'))
    fasta = confined_file(settings['fasta'], [references])
    files = {'reference.fa': fasta, 'reference.fa.fai': confined_file(str(fasta)+'.fai', [references])}
    with (output / 'samples.csv').open(encoding='utf-8', newline='') as handle:
        samples = list(csv.DictReader(handle))
    for i, sample in enumerate(samples):
        bam = confined_file(sample['bam_path'], [data])
        index = Path(str(bam)+'.bai')
        if not index.is_file():
            index = bam.with_suffix('.bai')
        files[f'sample-{i}.bam'] = bam
        files[f'sample-{i}.bam.bai'] = confined_file(index, [data])
    return files, samples, settings


def annotations(output):
    source = json.loads((Path(output)/'gene_structure.json').read_text(encoding='utf-8'))
    lines = []
    for tx in source['transcripts']:
        exons = sorted(tx['exons'], key=lambda e:e['start_1based'])
        start, end = exons[0]['start_1based']-1, exons[-1]['end_1based']
        cds = [c for e in exons for c in e['cds']]
        thick_start = min((c['start_1based']-1 for c in cds), default=start)
        thick_end = max((c['end_1based'] for c in cds), default=start)
        lines.append('\t'.join(map(str,[tx['chrom'],start,end,tx['id'],0,tx['strand'],
            thick_start,thick_end,0,len(exons),
            ','.join(str(e['end_1based']-e['start_1based']+1) for e in exons),
            ','.join(str(e['start_1based']-1-start) for e in exons)])))
    return source, '\n'.join(lines)+'\n'


def configuration(output, jid, data, references):
    _, samples, settings = assets(output, data, references)
    source, _ = annotations(output)
    tx = source['transcripts'][0]
    start = max(1,min(e['start_1based'] for t in source['transcripts'] for e in t['exons'])-500)
    end = max(e['end_1based'] for t in source['transcripts'] for e in t['exons'])+500
    base = f'/api/jobs/{jid}/igv/'
    tracks = [dict(name='解析に使用した遺伝子注釈', type='annotation', format='bed',
                   url=base+'annotations.bed', indexed=False, displayMode='EXPANDED', order=0)]
    for i, s in enumerate(samples):
        tracks.append(dict(name=f"{s['metadata_label']} · {s['srx_id']}", type='alignment', format='bam',
            url=base+f'sample-{i}.bam', indexURL=base+f'sample-{i}.bam.bai',
            color='#256aa0' if s['metadata_label']=='WT' else '#a65b28',
            showCoverage=True, showAlignments=True, viewAsPairs=False,
            visibilityWindow=200000, height=240, order=i+1))
    return dict(reference=dict(id=settings['assembly'], name=settings['assembly'],
                fastaURL=base+'reference.fa', indexURL=base+'reference.fa.fai'),
                locus=f"{tx['chrom']}:{start}-{end}", tracks=tracks, showUCSCLink=False)


def byte_range(header, size):
    """Return an inclusive single byte range, or raise on malformed requests."""
    import re
    match = re.fullmatch(r'bytes=(\d*)-(\d*)', header or '')
    if not match or not any(match.groups()) or size == 0:
        raise ValueError('Invalid byte range')
    first, last = match.groups()
    if first:
        start, end = int(first), min(int(last), size-1) if last else size-1
    else:
        if int(last) == 0:
            raise ValueError('Invalid suffix range')
        start, end = max(0,size-int(last)),size-1
    if start >= size or start > end:
        raise ValueError('Unsatisfiable byte range')
    return start, end
