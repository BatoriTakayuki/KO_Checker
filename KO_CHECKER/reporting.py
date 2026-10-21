"""Conclusion-first RNA reports, including conservative exploratory statistics."""
from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path
import random
from statistics import mean


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) and value >= 0 else None
    except (TypeError, ValueError):
        return None


def permutation_test(wt, ko):
    """Two-sided |mean difference| test; labels must be exchangeable.

    Enumerate up to 50,000 partitions, otherwise sample 19,999 with seed 0.
    Randomized p values include the observed partition (plus-one correction).
    """
    wt, ko = sorted(wt), sorted(ko)
    values = wt + ko
    n_wt, n_ko = len(wt), len(ko)
    total = math.fsum(values)
    observed = abs(mean(ko) - mean(wt))
    tolerance = max(1.0, max(values)) * 1e-12
    partitions = math.comb(len(values), n_ko)
    exact = partitions <= 50000
    trials = partitions if exact else 19999
    rng = random.Random(0)
    choices = (itertools.combinations(range(len(values)), n_ko) if exact else
               (rng.sample(range(len(values)), n_ko) for _ in range(trials)))
    extreme = 0
    for indices in choices:
        selected = math.fsum(values[i] for i in indices)
        statistic = abs(selected / n_ko - (total - selected) / n_wt)
        extreme += statistic >= observed - tolerance
    p = extreme / trials if exact else (extreme + 1) / (trials + 1)
    return {'status': 'tested', 'method': 'exact' if exact else 'monte_carlo',
            'p_value': p, 'alpha': 0.05, 'significant': p < 0.05,
            'partitions': trials, 'seed': None if exact else 0}


def expression_summary(samples, confirmed):
    groups = {'WT': [], 'KO': []}
    excluded = []
    for row in samples:
        label = row.get('metadata_label')
        if label not in groups:
            continue
        cpm = number(row.get('gene_unique_fragment_cpm'))
        library = number(row.get('unique_primary_library_fragments'))
        if cpm is None or library is None or library <= 0:
            excluded.append(row.get('srx_id', ''))
        else:
            groups[label].append((row.get('srx_id', ''), cpm))
    wt, ko = ([v for _, v in groups[label]] for label in ('WT', 'KO'))
    wt_mean, ko_mean = mean(wt) if wt else None, mean(ko) if ko else None
    ratio = ko_mean / wt_mean if wt_mean and ko_mean is not None else None
    result = {'status': 'compared' if ratio is not None else 'unavailable',
              'wt_n': len(wt), 'ko_n': len(ko), 'wt_mean_cpm': wt_mean,
              'ko_mean_cpm': ko_mean, 'ratio': ratio,
              'decrease_percent': 100 * (1 - ratio) if ratio is not None else None,
              'excluded_samples': excluded,
              'baseline': 'all_WT_labelled_samples_with_valid_library',
              'samples': [{'srx_id': sid, 'metadata_label': label, 'cpm': cpm,
                           'ratio': cpm / wt_mean if wt_mean else None}
                          for label, rows in groups.items() for sid, cpm in rows]}
    if ratio is None:
        text = ('WT群の発現量が0のため、WTに対する低下率を計算できません。' if wt_mean == 0 else
                'WT・KO両群の有効なライブラリ計数が揃わず、発現量を比較できません。')
    else:
        change = abs(result['decrease_percent'])
        direction = '低下' if ratio < 1 else '増加' if ratio > 1 else '変化なし'
        text = (f'KO群の相対発現量（CPM）はWT群の{ratio * 100:.1f}%です'
                f'（{change:.1f}%{direction}、WT n={len(wt)} / KO n={len(ko)}）。')
    if excluded:
        text += ' ライブラリ計数が無効な検体を比較から除外しました: ' + ', '.join(excluded) + '。'
    result['text'] = text
    if ratio is None:
        test = {'status': 'unavailable', 'text': '有意差は評価できません。'}
    elif min(len(wt), len(ko)) < 2:
        test = {'status': 'insufficient_replicates', 'text': '各群2検体以上の独立した生物学的反復がないため、有意差は評価できません。'}
    elif not confirmed:
        test = {'status': 'not_confirmed', 'text': '有意差は未評価です。独立した生物学的反復と、群間で比較可能な実験条件の確認が必要です。'}
    else:
        test = permutation_test(wt, ko)
        verdict = '有意差あり' if test['significant'] else '有意差を検出せず（差がないことの証明ではありません）'
        mode = '全組合せ' if test['method'] == 'exact' else '乱数近似'
        test['text'] = (f'探索的なCPM平均差の両側置換検定: {verdict}'
                        f'（p={test["p_value"]:.4g}、{mode} {test["partitions"]:,}通り、未補正p、基準0.05）。')
    result['test'] = test
    return result


def build_conclusion(tables, confirmed_replicates=False, thresholds=None):
    thresholds = thresholds or {}
    samples = tables.get('samples', [])
    ko = [s for s in samples if s.get('metadata_label') == 'KO']
    exons = tables.get('exon_coverage', [])
    ko_exons = [r for r in exons if r.get('metadata_label') == 'KO']
    region_models = tables.get('region_deletions', [])
    region_ids = {r['srx_id'] for r in region_models if r.get('ko_support') == 'yes'}
    counts = {state: sum(s.get('predicted_class') == state for s in ko)
              for state in ('KO_like', 'WT_like', 'undetermined')}
    low = {s['srx_id'] for s in ko if s.get('coverage_warning')}
    low.update(r['srx_id'] for r in ko_exons
               if number(r.get('sufficient_depth_fraction')) is not None
               and float(r['sufficient_depth_fraction']) < .8)
    if counts['KO_like']:
        state = 'ko_supported'
        text = f'KO群{len(ko)}検体中{counts["KO_like"]}検体で、KOを支持するRNA上の配列・スプライス異常またはコード領域欠損が検出されました。'
        if region_ids:
            text += ' ' + ', '.join(sorted(region_ids)) + 'はKO様（RNA領域欠損）です。欠損後の配列は注釈エクソン境界を仮定した予測で、DNA欠失は未確定です。'
    elif low:
        state = 'insufficient_coverage'
        text = 'KOを支持する配列異常を示せませんでした。' + ', '.join(sorted(low)) + 'では読み取りが不足する領域があり、変異の有無を判断できません。'
    elif counts['WT_like']:
        state = 'reference_supported'
        text = f'KO群{len(ko)}検体中{counts["WT_like"]}検体で、観測したKO候補位置の参照配列または正常接合を支持します。'
    elif ko:
        state = 'undetermined'
        text = '観測範囲ではKOを支持する配列異常を検出できませんでした。配列からの結論は判定保留です。'
    else:
        state = 'unavailable'
        text = 'KO群の検体がなく、KO群の配列評価をまとめられません。'
    if counts['KO_like'] and (counts['undetermined'] or counts['WT_like']):
        text += f' その他はWT様の局所証拠{counts["WT_like"]}検体、判定保留{counts["undetermined"]}検体です。'
    sequence = {'status': state, 'text': text, 'counts': counts,
                'low_coverage_samples': sorted(low),
                'note': 'RNAで変異が見えないことはKO不成立を意味しません。DNA上のKO成立やタンパク質消失、遺伝子全体のWTは確定しません。'}
    expression = expression_summary(samples, confirmed_replicates)
    depleted = [r for r in ko_exons if r.get('comparison_status') == 'suspected_local_depletion']
    insufficient = sum(r.get('comparison_status', '').startswith('insufficient_') for r in ko_exons)
    strongest = sorted(depleted, key=lambda r: number(r.get('relative_depth_vs_wt')) or 0)[:10]
    if depleted:
        exon_text = (f'遺伝子内の他のエクソンで補正した、局所的な読み取り低下の所見が'
                     f'{len(depleted)}件（検体×エクソン）ありました。統計検定による判定ではありません。')
    elif ko_exons:
        exon_text = '評価できた範囲で、局所的な読み取り低下の基準を満たすエクソンはありませんでした。'
    else:
        exon_text = 'エクソン比較データがありません。Webで再解析するとエクソンごとの比較を追加できます。'
    if insufficient:
        exon_text += f' 深度不足などにより評価できない所見は{insufficient}件です。'
    if region_models:
        exon_text = (f'WTとのライブラリ補正比較と残存領域の対照に基づくRNA領域欠損モデルが{len(region_models)}件（検体×転写産物）あります。'
                     f' 従来の他エクソン補正による比較では、深度条件などで{insufficient}件が保留ですが、新しい領域欠損評価は別に実施しています。')
    exon_summary = {'text': exon_text, 'depleted_count': len(depleted),
                    'insufficient_count': insufficient, 'strongest': strongest}
    evidence = {(r.get('srx_id'), r.get('event_id'))
                for name in ('events', 'junctions') for r in tables.get(name, [])
                if r.get('comparison_status') == 'enriched_vs_comparators'}
    ko_ids = {s.get('srx_id') for s in ko}
    nmd_samples = {r.get('srx_id') for r in tables.get('consequences', [])
                   if r.get('srx_id') in ko_ids and r.get('nmd_prediction') == 'possible'
                   and r.get('effect') in ('frameshift', 'stop_gained', 'exon_skip_frameshift')
                   and (r.get('srx_id'), r.get('event_id')) in evidence}
    lower = {s['srx_id'] for s in expression['samples']
             if s['metadata_label'] == 'KO' and s['ratio'] is not None and s['ratio'] < 1}
    associated = sorted(nmd_samples & lower)
    reduced = expression['ratio'] is not None and expression['ratio'] < 1
    if associated and reduced:
        nmd_state = 'consistent'
        nmd_text = ('NMDが起こり得る配列異常の予測と発現低下が同じ検体にあります（'
                    + ', '.join(associated) + '）。NMDの関与と整合しますが、機序は未確認です。')
    elif reduced:
        nmd_state = 'unresolved'
        nmd_text = '発現低下を認めます。NMDを含むRNA分解の促進、転写低下などの原因は、このデータだけでは区別できません。'
    elif nmd_samples:
        nmd_state = 'sequence_prediction_only'
        nmd_text = 'NMDが起こり得る配列上の予測がありますが、同じ検体での発現低下と群全体の低下は揃っていません。'
    else:
        nmd_state = 'unavailable'
        nmd_text = 'NMDの関与を示唆する配列予測と発現低下の組合せは得られていません。'
    if region_ids:
        predictions = [r.get('nmd_prediction', 'unknown') for r in region_models if r.get('ko_support') == 'yes']
        possible = predictions.count('possible_under_model')
        unlikely = predictions.count('unlikely_under_model')
        nmd_text += (f' 領域欠損を仮定した主要CDSモデルでは、NMDが起こり得る予測{possible}件、'
                     f'起こりにくい予測{unlikely}件、その他・不明{len(predictions)-possible-unlikely}件です。'
                     '終止コドンと残存接合の位置による条件付き予測で、NMDの実測ではありません。')
    overall = ('配列による判定が保留でも、対象遺伝子の相対発現低下は独立した所見です。タンパク質量や機能の低下の程度は未確認です。'
               if reduced else '配列・発現量・エクソンの所見を合わせて確認してください。RNA解析による暫定評価です。')
    if region_ids:
        overall = '主要なコード配列の領域欠損を支持するRNA所見からKO様と推定します。全体のCPMが保たれていても、正常な全長mRNAの維持を意味しません。実際の接合境界・タンパク質とDNA欠失は未確認です。'
    return {'version': 1, 'gene': next((s.get('gene') for s in samples if s.get('gene')), ''),
            'sequence': sequence, 'expression': expression, 'exons': exon_summary,
            'nmd': {'status': nmd_state, 'text': nmd_text, 'samples': associated},
            'overall': overall,
            'region_loss': {'model_count': len(region_models), 'ko_samples': sorted(region_ids),
                            'text': 'RNA領域欠損の条件付き転写産物・アミノ酸配列は領域欠損表とFASTAを参照してください。' if region_models else 'RNA領域欠損モデルはありません。'},
            'notes': ['発現量はライブラリ断片数で補正した相対指標CPMで、細胞あたりの絶対mRNA量ではありません。',
                      '置換検定は独立・非対応で群ラベルを交換可能とする探索的検定です。バッチ補正やRNA-seqの組成補正、全遺伝子の発現差解析は行いません。',
                      'gRNA標的位置は未指定です。カバレッジは注釈エクソンの評価で、標的塩基の観測を保証しません。'],
            'min_event_depth': thresholds.get('min_event_depth', 10)}


def read_conclusion(output, confirmed_replicates=False):
    """Read full CSVs, so the 1,000-row browser limit never truncates conclusions."""
    output = Path(output)
    tables = {}
    for name in ('samples', 'events', 'junctions', 'consequences', 'exon_coverage', 'region_deletions'):
        path = output / (name + '.csv')
        tables[name] = []
        if path.is_file():
            with path.open(newline='', encoding='utf-8') as handle:
                tables[name] = [{k: v for k, v in row.items()
                                 if k not in ('reference_protein', 'predicted_protein', 'predicted_mrna')}
                                for row in csv.DictReader(handle)]
    try:
        settings = json.loads((output / 'run.json').read_text(encoding='utf-8'))
    except (OSError, ValueError):
        settings = {}
    return build_conclusion(tables, confirmed_replicates, settings.get('thresholds'))


def conclusion_text(report):
    expression = report['expression']
    highlights = '\n'.join(
        f'  {r["srx_id"]}: {r["chrom"]}:{r["exon_start_1based"]}-{r["exon_end_1based"]}'
        f'（他エクソン補正後のWT比 {float(r["relative_depth_vs_wt"]) * 100:.1f}%）'
        for r in report['exons']['strongest'])
    return '\n\n'.join([
        (report.get('gene') or '対象遺伝子') + ' — 解析の結論',
        '配列による評価: ' + report['sequence']['text'],
        '発現量: ' + expression['text'] + ' ' + expression['test']['text'],
        'エクソン: ' + report['exons']['text'] + ('\n' + highlights if highlights else ''),
        'NMDの解釈: ' + report['nmd']['text'],
        '総合評価: ' + report['overall'], report['sequence']['note'], *report['notes'],
        report['region_loss']['text'],
    ]) + '\n'


def write_conclusion(output, confirmed_replicates=False):
    report = read_conclusion(output, confirmed_replicates)
    for name, content in [('conclusion.json', json.dumps(report, ensure_ascii=False, indent=2) + '\n'),
                          ('conclusion.txt', conclusion_text(report))]:
        target = Path(output) / name
        temporary = target.with_suffix(target.suffix + '.partial')
        temporary.write_text(content, encoding='utf-8')
        temporary.replace(target)
    return report
