"""Reference-coordinate transcript diagrams with observed and conditional evidence."""
import csv
import json
from pathlib import Path

import rnako


def table(output, name):
    path = Path(output) / (name + '.csv')
    if not path.is_file():
        return []
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def codon(sequence, positions, offset, status):
    if offset is None or offset < 0 or offset+3 > len(sequence):
        return None
    return {'bases': sequence[offset:offset+3], 'cdna_1based': offset+1,
            'genomic_positions_1based': [p+1 for p in positions[offset:offset+3]], 'status': status}


def coordinate_template(transcript, fa):
    seq, mapping, start, end, _ = rnako.template(transcript, fa, rnako.resolve_reference_name(transcript.chrom, fa))
    positions = [p for p, index in sorted(mapping.items(), key=lambda item: item[1])]
    return seq, positions, start, end


def coding_segments(positions, start, length):
    """Project translated bases to the reference, including splice boundaries."""
    if start is None or length is None:
        return None
    runs = []
    for position in positions[start:start+length*3]:
        if runs and abs(position-runs[-1][-1]) <= 1:
            runs[-1].append(position)
        else:
            runs.append([position])
    return [{'start_1based':min(run)+1,'end_1based':max(run)+1} for run in runs]


def reference_model(transcript, fa):
    seq, positions, start, end = coordinate_template(transcript, fa)
    valid_start = start is not None and seq[start:start+3] == 'ATG'
    peptide, stop = rnako.translate(seq[start:]) if valid_start else ('', None)
    return {'kind': 'reference_annotation', 'exon_numbers': list(range(1, len(transcript.exons)+1)),
            'start': codon(seq, positions, start, 'reference_annotation' if valid_start else 'noncanonical_or_partial_annotation'),
            'stop': codon(seq, positions, start+stop if stop is not None else
                          end if end is not None and seq[end:end+3] in ('TAA','TGA','TAG') else None,
                          'reference_annotation'),
            'transcript_length': len(seq), 'protein_length': len(peptide) if valid_start else None}


def loss_model(transcript, row, fa):
    seq, positions, start, end = coordinate_template(transcript, fa)
    removed = {int(x) for x in row['lost_exons'].split(';') if x}
    intervals = [ex for i, ex in enumerate(transcript.ordered_exons(), 1) if i in removed]
    keep = [i for i, p in enumerate(positions) if not any(a <= p < b for a,b in intervals)]
    new_seq = ''.join(seq[i] for i in keep)
    new_positions = [positions[i] for i in keep]
    kept_start = start is not None and all(i in keep for i in range(start,start+3))
    new_start = keep.index(start) if kept_start else None
    stop = int(row['stop_codon_cdna_1based'])-1 if row.get('stop_codon_cdna_1based') else None
    return {'kind': 'region_loss', 'event_id': row['event_id'], 'evidence': row['evidence'],
            'effect': row['effect'], 'removed_exon_numbers': sorted(removed),
            'exon_numbers': [i for i in range(1,len(transcript.exons)+1) if i not in removed],
            'start': codon(new_seq,new_positions,new_start,'conditional_model') if kept_start else None,
            'stop': codon(new_seq,new_positions,stop,'conditional_model'),
            'transcript_length': len(new_seq),
            'protein_length': len(row['predicted_protein']) if row.get('predicted_protein') else None,
            'coding_segments': coding_segments(new_positions,new_start,len(row['predicted_protein'])) if row.get('predicted_protein') else None,
            'prediction_status': row['protein_prediction_status'], 'ko_support': row['ko_support'],
            'note': row['assumptions']}


def variant_model(transcript, event, fa):
    seq, positions, start, end = coordinate_template(transcript, fa)
    mapping = {p:i for i,p in enumerate(positions)}
    mutation = rnako.altered_transcript(seq,mapping,transcript,
        {'pos':int(event['pos_1based']),'ref':event['ref'],'alt':event['alt']})
    if mutation is None or start is None:
        return None
    new_seq, edit, removed, added = mutation
    # Inserted bases have no reference coordinate; place their marker at the
    # insertion anchor without pretending they occupy new genomic positions.
    new_positions = positions[:edit] + [int(event['pos_1based'])-1]*added + positions[edit+removed:]
    effect = rnako.protein_effect(transcript,
        {'pos':int(event['pos_1based']),'ref':event['ref'],'alt':event['alt']},fa,
        rnako.resolve_reference_name(transcript.chrom,fa))
    start_affected = edit < start+3 and edit+max(removed,1) > start
    new_start = start + (added-removed if edit+removed <= start else 0)
    valid_start = not start_affected and new_seq[new_start:new_start+3] == 'ATG'
    protein, stop = rnako.translate(new_seq[new_start:]) if valid_start else ('',None)
    return {'kind':'single_variant', 'event_id':event['event_id'], 'effect':effect['effect'],
            'exon_numbers': list(range(1,len(transcript.exons)+1)), 'removed_exon_numbers':[],
            'start':codon(new_seq,new_positions,new_start,'start_changed' if start_affected else 'single_variant_model'),
            'stop':codon(new_seq,new_positions,new_start+stop if stop is not None else None,'single_variant_model'),
            'protein_length':len(protein) if valid_start else None, 'transcript_length':len(new_seq),
            'coding_segments':coding_segments(new_positions,new_start,len(protein)) if valid_start else None,
            'note':'1変異のみを適用した条件付きモデル。複数変異の位相や検体全体の配列を確定しません。'}


def write_structure(output, fasta, gtf, gene, transcripts=None):
    output = Path(output)
    transcripts = transcripts if transcripts is not None else rnako.load_gene(gtf,gene)
    samples = table(output,'samples')
    observations = table(output,'observations')
    events = table(output,'events')
    # CLI does not write observations.csv; events still carry local observations.
    if not observations:
        observations = events
    coverage = table(output,'exon_coverage')
    models = table(output,'region_deletions')
    structures = []
    with rnako.pysam.FastaFile(str(fasta)) as fa:
        for t in transcripts.values():
            exons = [{'number':i,'start_1based':a+1,'end_1based':b,
                      'cds':[{'start_1based':max(a,c)+1,'end_1based':min(b,d)}
                             for c,d in t.cds if max(a,c)<min(b,d)]}
                     for i,(a,b) in enumerate(t.ordered_exons(),1)]
            loss = {r['srx_id']:loss_model(t,r,fa) for r in models if r['transcript_id']==t.tid}
            variants = []
            for event in events:
                if any(a<=int(event['pos_1based'])-1<b for a,b in t.exons):
                    prediction = variant_model(t,event,fa)
                    if prediction:
                        variants.append({'srx_id':event['srx_id'],**prediction})
            structures.append({'id':t.tid,'chrom':t.chrom,'strand':t.strand,'exons':exons,
                'cds_length':sum(b-a for a,b in t.cds), 'reference':reference_model(t,fa),
                'loss_models':loss,'variant_models':variants})
    structures.sort(key=lambda t:(-any(m['ko_support']=='yes' for m in t['loss_models'].values()),
                                  -(t['reference']['start'] is not None and t['reference']['start']['bases']=='ATG'),
                                  -t['cds_length'],t['id']))
    result = {'version':1,'gene':gene,'coordinate_system':'1-based inclusive; transcript orientation 5prime to 3prime',
              'samples':[{'id':s['srx_id'],'group':s['metadata_label']} for s in samples],
              'transcripts':structures,'coverage':coverage,'variants':observations,
              'note':'WTは参照注釈に観測値を重ねた図。KOの欠損・変異モデルは条件付き予測。イントロンを短縮した模式図はゲノム距離に比例しません。'}
    target = output / 'gene_structure.json'
    temporary = target.with_suffix('.json.partial')
    temporary.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    temporary.replace(target)
    return result
