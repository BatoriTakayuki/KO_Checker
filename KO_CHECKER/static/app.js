'use strict';
const exonStatus = {no_coverage:'読み取りなし',low_coverage:'深度不足',covered:'深度基準を満たす',reference_sample:'WT比較用検体',insufficient_wt_coverage:'WTの深度不足で評価不能',insufficient_sample_coverage:'検体の深度不足で評価不能',suspected_local_depletion:'局所的な読み取り低下の疑い',no_specific_depletion:'局所低下の基準に該当せず'};
const $ = id => document.getElementById(id);
const statusText = {queued:'待機中',running:'解析中',paused:'中断中',completed:'完了',failed:'失敗'};
const stages = {queued:'実行待ち',preflight:'入力確認',reference:'参照・索引',reference_download:'参照ダウンロード',reference_extract:'参照の展開',reference_index:'FASTA索引',reference_annotation:'遺伝子注釈の確認',hisat2_index:'HISAT2索引の構築',download:'SRA / FASTQ',sra_download:'SRAの取得',sra_validate:'SRAの整合性検査',fastq_convert:'FASTQへの変換',mapping:'HISAT2',merge:'SRX単位の統合',analysis:'RNA評価',completed:'完了',failed:'失敗',interrupted:'中断'};
const explanations = {
  'Insufficient uniquely mapped exonic fragments to interpret lack of evidence':'対象遺伝子の有効断片が少なく、判断できません。',
  'Potential damaging event enriched relative to WT-labelled samples; RNA evidence only':'比較検体より多い有害候補変異、または未注釈エクソンスキップを検出しました（RNA上の証拠）。',
  'No KO-specific disruptive event demonstrated; cannot establish wild-type genotype from RNA':'KOを支持する特異的な異常を示せません。RNAだけではWTを確定できません。',
  'Reference bases at covered KO-like variant site(s) and/or annotated junction support; other lesions remain possible':'確認できたKO候補位置の変異率が低い、または正常なスプライス接合を支持します。他の異常は否定できません。'
};
let active = null, rendered = null, currentJob = null, resultData = null, polling = false;
async function api(path, options) {
  const res = await fetch(path, options);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
  return data;
}
function el(tag, text, className) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (className) node.className = className;
  return node;
}
function number(value, digits=3) {
  return value === '' || value === null || value === undefined ? '—' : Number(value).toLocaleString('ja-JP',{maximumFractionDigits:digits});
}
function choose(id) {
  active = id; rendered = null; currentJob = null; resultData = null;
  location.hash = id; $('results').hidden = true; $('empty').hidden = true; $('job').hidden = false;
  refresh();
}
async function submit(payload) {
  $('submit').disabled = true; $('form-error').textContent = '';
  try {
    const job = await api('/api/jobs',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});
    choose(job.id);
  } catch (err) { $('form-error').textContent = err.message; }
  finally { $('submit').disabled = false; }
}
$('analysis-form').addEventListener('submit', event => {
  event.preventDefault(); const data=Object.fromEntries(new FormData(event.target));
  data.pause_after_reference=$('pause-after-reference').checked;submit(data);
});
async function control(operation) {
  if (!active) return;
  const selected=active;$('pause').disabled=true;$('resume').disabled=true;
  try {
    await api(`/api/jobs/${selected}/${operation}`,{method:'POST',headers:{'Content-Type':'application/json'},body:'{}'});
    await refresh();
  } catch(err) {$('form-error').textContent=err.message;}
  finally {$('pause').disabled=false;$('resume').disabled=false;}
}
$('pause').addEventListener('click',()=>control('pause'));
$('resume').addEventListener('click',()=>control('resume'));
function duration(seconds) {
  seconds=Math.max(0,Math.floor(seconds));
  return `${Math.floor(seconds/3600)}時間${Math.floor(seconds%3600/60)}分${seconds%60}秒`;
}
function progressDetails(job) {
  const start=Date.parse(job.attempt_started_at || job.created_at);
  const end=job.status==='running'?Date.now():Date.parse(job.paused_at || job.updated_at);
  $('elapsed').textContent=`今回の実行時間: ${duration((end-start)/1000)}`;
  $('run-count').textContent=job.runs_total?`run完了: ${job.runs_done || 0} / ${job.runs_total}`:'';
  $('last-update').textContent=`状態更新: ${new Date(job.updated_at).toLocaleTimeString('ja-JP')}`;
  $('activity').textContent=job.activity || (job.active_command?`${job.active_command}の終了を待っています。`:'');
  $('storage-progress').textContent=job.storage_bytes!==null && job.storage_bytes!==undefined?`${job.storage_label}: ${number(job.storage_bytes/1000000,1)} MB`:'';
  $('pause').hidden=!['queued','running'].includes(job.status);
  $('pause').disabled=!!job.pause_requested;
  $('resume').hidden=job.status!=='paused';
  $('pause-info').textContent=job.pause_requested?'中断中の表示になるまで、DockerとPCを終了しないでください。通信待ちなどで数十秒かかる場合があります。':job.status==='paused'?'中断済みです。完了済みの参照・FASTQ・BAMを保持しています。再開時は未完了の計算工程だけをやり直します。':['queued','running'].includes(job.status)?'中断すると実行中の処理を終了します。参照取得は配信元が対応していれば続きから再開し、計算途中の索引・マッピングは工程の最初からやり直します。':'';
  const transfer=job.transfer;$('transfer').hidden=!transfer;
  if(transfer) {
    const total=Number(transfer.total_bytes),done=Number(transfer.done_bytes),rate=Number(transfer.rate_bps);
    if(total>0)$('transfer-progress').value=Math.min(100,done/total*100);else $('transfer-progress').removeAttribute('value');
    const size=`${number(done/1000000,1)} MB`+(total>0?` / ${number(total/1000000,1)} MB (${number(done/total*100,1)}%)`:'');
    const speed=rate>0&&job.status==='running'?` · ${number(rate/1000000,2)} MB/s`:'';
    const eta=rate>0&&total>done&&job.status==='running'?` · このファイルの残り目安 ${duration((total-done)/rate)}`:'';
    $('transfer-text').textContent=`${transfer.name}: ${size}${speed}${eta}`+(transfer.resumed_from>0?` · ${number(transfer.resumed_from/1000000,1)} MBから再開`:'');
  }
}
$('retry').addEventListener('click', () => {
  if (currentJob) submit({...currentJob.request,wt:currentJob.request.wt.join('\n'),ko:currentJob.request.ko.join('\n')});
});
function table(id, rows, columns) {
  const target = $(id); target.replaceChildren();
  const head = el('tr');
  columns.forEach(([,label]) => head.append(el('th',label)));
  const thead = el('thead'); thead.append(head); target.append(thead);
  const body = el('tbody');
  for (const row of rows) {
    const tr = el('tr');
    for (const [key] of columns) {
      const td = el('td'); const value = row[key] ?? '';
      if (key === 'predicted_class') td.append(el('span',value,`badge ${value}`));
      else if (key === 'transcript_exons' && value) {
        const detail=el('details');const entries=value.split(';');detail.append(el('summary',`転写産物 ${entries.length}件`),el('p',entries.join('\n')));td.append(detail);
      }
      else if (key.endsWith('_protein') && value) {
        const detail = el('details'); detail.append(el('summary',`${value.length} aa · 配列を表示`),el('pre',value)); td.append(detail);
      } else td.textContent = value === '' ? '—' : (key === 'classification_reason' ? (explanations[value] || value) : ((key === 'coverage_status' || key === 'comparison_status') ? (exonStatus[value] || value) : String(value)));
      tr.append(td);
    }
    body.append(tr);
  }
  if (!rows.length) {const tr=el('tr'); const td=el('td','該当するデータはありません。');td.colSpan=columns.length;tr.append(td);body.append(tr);}
  target.append(body);
}
function bars(id, rows, maximum, fraction=false) {
  const target = $(id); target.replaceChildren();
  for (const row of rows) {
    const line=el('div',undefined,'bar-row');
    const label=el('span',row.id,'bar-label');label.title=`${row.id} (${row.group})`;
    const track=el('div',undefined,'bar-track');
    const bar=el('div',undefined,`bar ${row.group === 'WT'?'wt':'ko'}`);
    const valid=row.value !== '' && row.value !== undefined && row.value !== null && Number.isFinite(Number(row.value));
    bar.style.width = valid ? `${Math.min(100,Math.max(0,Number(row.value)/(maximum||1)*100))}%` : '0%';track.append(bar);
    line.append(label,track,el('span',valid ? number(Number(row.value)*(fraction?100:1),2)+(fraction?'%':'') : '未掲載','bar-value'));
    target.append(line);
  }
}
function eventChart() {
  if (!resultData) return;
  const key=$('event-select').value;
  bars('event-chart',resultData.tables.samples.map(s=>({id:s.srx_id,group:s.metadata_label,
    value:resultData.tables.observations.find(e=>e.srx_id===s.srx_id && `${e.chrom}:${e.pos_1based}:${e.ref}:${e.alt}`===key)?.local_alt_fraction})),1,true);
  $('event-depth').textContent=resultData.tables.observations.filter(e=>`${e.chrom}:${e.pos_1based}:${e.ref}:${e.alt}`===key).map(e=>`${e.srx_id}: ${e.alt_reads}/${e.depth_at_anchor} ${e.coverage_status==='low_depth'?'(低深度)':''}`).join(' / ');
}
$('event-select').addEventListener('change',eventChart);
function results(data, job) {
  resultData=data;const t=data.tables;
  $('results').hidden=false;
  $('summary').replaceChildren();
  for (const state of ['KO_like','WT_like','undetermined']) {
    const box=el('div',state);box.prepend(el('strong',t.samples.filter(s=>s.predicted_class===state).length));$('summary').append(box);
  }
  table('samples-table',t.samples,[['srx_id','入力ID'],['metadata_label','指定群'],['predicted_class','RNA評価'],['gene_unique_fragment_cpm','CPM'],['expression_relative_to_wt','WT平均比'],['coverage_warning','カバレッジ注意'],['classification_reason','判定根拠']]);
  const exons=t.exon_coverage || [];
  const percent=value=>`${number(Number(value)*100,1)}%`;
  const exonRows=exons.map(row=>({...row,region:`${row.chrom}:${row.exon_start_1based}-${row.exon_end_1based}`,covered_percent:percent(row.covered_fraction),sufficient_percent:percent(row.sufficient_depth_fraction)}));
  table('exons-table',exonRows,[['srx_id','入力ID'],['metadata_label','指定群'],['region','エクソン座標 (1-based)'],['coverage_status','読み取り状態'],['comparison_status','WT比較の評価'],['mean_depth','平均深度'],['overlapping_fragments','重なる断片数'],['min_depth','最低深度'],['covered_percent','読まれた塩基'],['sufficient_percent','深度10以上の塩基'],['relative_depth_vs_wt','他エクソン補正後のWT比'],['wt_reference_samples','比較に使用したWT数'],['skipping_junction_fragments','スキップ接合支持断片'],['longest_zero_depth_region','最長ゼロ深度区間'],['transcript_exons','転写産物・エクソン番号']]);
  const koExons=exons.filter(row=>row.metadata_label==='KO');
  const suspected=koExons.filter(row=>row.comparison_status==='suspected_local_depletion').length;
  const insufficient=koExons.filter(row=>row.comparison_status.startsWith('insufficient_')).length;
  $('exon-summary').textContent=exons.length?`KO群の検体×エクソン ${koExons.length}件中、局所的な読み取り低下の疑い ${suspected}件、深度不足で評価不能 ${insufficient}件。RNA評価（KO_like等）とは別の補助所見です。`:'この解析にはエクソン比較データがありません。「同じ条件で再実行」で追加できます。';
  const cpm=t.samples.map(s=>({id:s.srx_id,group:s.metadata_label,value:s.gene_unique_fragment_cpm}));
  bars('cpm-chart',cpm,Math.max(1,...cpm.map(r=>Number(r.value))));
  const ratios=t.samples.map(s=>({id:s.srx_id,group:s.metadata_label,value:s.expression_relative_to_wt}));
  bars('ratio-chart',ratios,Math.max(1,...ratios.map(r=>Number(r.value))));
  // Include a direct WT/KO mean comparison alongside individual samples.
  const mean=el('p',undefined,'hint');
  mean.textContent=['WT','KO'].map(group=>{const values=t.samples.filter(s=>s.metadata_label===group && s.gene_unique_fragment_cpm!=='').map(s=>Number(s.gene_unique_fragment_cpm));return `${group}平均CPM: ${values.length?number(values.reduce((a,b)=>a+b,0)/values.length):'—'} (n=${values.length})`;}).join(' / ');
  $('cpm-chart').append(mean);
  const options=new Map();
  t.observations.forEach(e=>options.set(`${e.chrom}:${e.pos_1based}:${e.ref}:${e.alt}`,`${e.chrom}:${e.pos_1based} ${e.ref} → ${e.alt}`));
  $('event-select').replaceChildren();
  for (const [value,label] of options) {const opt=el('option',label);opt.value=value;$('event-select').append(opt);}
  if (!options.size) $('event-select').append(el('option','変異候補は検出されませんでした'));
  eventChart();
  table('events-table',t.events,[['srx_id','入力ID'],['chrom','染色体'],['pos_1based','位置 (1-based)'],['ref','参照'],['alt','変異'],['alt_reads','支持リード'],['depth_at_anchor','深度'],['local_alt_fraction','局所変異率'],['wt_label_max_fraction','比較側の最大率'],['comparison_status','比較結果']]);
  table('junctions-table',t.junctions,[['srx_id','入力ID'],['intron_start_1based','開始'],['intron_end_1based','終了'],['support_fragments','支持断片'],['annotation_status','注釈'],['local_junction_depth','局所分母'],['local_skip_fraction','局所支持率'],['comparator_max_skip_fraction','比較側の最大率'],['comparison_status','比較結果'],['skipped_transcripts','スキップ対象']]);
  table('effects-table',t.consequences,[['srx_id','入力ID'],['event_id','イベントID'],['transcript_id','転写産物'],['effect','影響'],['first_changed_aa','最初の変化 (aa)'],['reference_protein','参照タンパク質'],['predicted_protein','予測タンパク質'],['stop_to_last_junction_nt','終止→最終接合 (nt)'],['nmd_prediction','NMD予測']]);
  $('table-limit').textContent=Object.entries(data.counts).map(([name,n])=>`${name}: ${n}件`).join(' / ')+'。各表は先頭1,000件まで表示します。全データはCSVから取得できます。';
  $('downloads').replaceChildren();
  for (const name of ['samples.csv','events.csv','junctions.csv','consequences.csv','observations.csv',...(exons.length?['exon_coverage.csv']:[]),'run.json','provenance.json','input_samples.csv','pipeline.log']) {
    const link=el('a',name);link.href=`/api/jobs/${job.id}/files/${name}`;link.download=name;$('downloads').append(link);
  }
}
async function refresh() {
  if (polling) return;polling=true;
  try {
    const jobs=await api('/api/jobs');$('history').replaceChildren();
    for (const job of jobs) {
      const button=el('button',`${job.request.gene} · ${job.request.assembly} / ${statusText[job.status]}`);
      button.append(el('small',new Date(job.created_at).toLocaleString('ja-JP')));
      button.addEventListener('click',()=>choose(job.id));$('history').append(button);
    }
    if (!jobs.length) $('history').append(el('p','まだ解析はありません。','muted'));
    if (!active && jobs.length) {active=jobs[0].id;$('empty').hidden=true;$('job').hidden=false;}
    if (!active) return;
    const selected=active;
    const job=await api(`/api/jobs/${selected}`);
    if (selected!==active) return;
    currentJob=job;
    $('job-id').textContent=`ANALYSIS · ${job.id.slice(0,8)}`;
    $('job-title').textContent=`${job.request.gene} / ${job.request.assembly}${job.threads ? ` / ${job.threads}スレッド` : ''}`;
    $('status').textContent=job.pause_requested?'中断処理中':statusText[job.status];$('status').className=`badge ${job.status}`;
    progressDetails(job);
    $('job-message').textContent=job.message;$('progress').value=job.progress;
    $('percent').textContent=`${job.progress}%（工程の目安）`;$('stage').textContent=stages[job.stage]||job.stage;
    $('job-error').textContent=job.status==='failed'?'処理が停止しました。ログで原因を確認し、解消後に再実行してください。':'';
    $('retry').hidden=job.status!=='failed';$('log-download').href=`/api/jobs/${job.id}/files/pipeline.log`;
    if ($('log-panel').open || job.status==='failed') {
      const log=await api(`/api/jobs/${selected}/log`);if(selected===active)$('log').textContent=log.text||'ログはまだありません。';
    }
    if (job.status==='completed' && rendered!==selected) {
      const data=await api(`/api/jobs/${selected}/results`);
      if(selected===active){results(data,job);rendered=selected;}
    }
  } catch(err) {$('form-error').textContent=`通信エラー: ${err.message}`;}
  finally {polling=false;}
}
$('log-panel').addEventListener('toggle',()=>{if($('log-panel').open)refresh();});
async function init() {
  try {
    const health=await api('/api/health');
    if (health.default_threads) $('threads').placeholder=`既定値: ${health.default_threads}スレッド（空欄可）`;
    const missing=Object.entries(health.tools).filter(([,ok])=>!ok).map(([name])=>name);if(!health.pysam)missing.push('pysam');
    $('health').textContent=missing.length?`解析環境が不足しています: ${missing.join(', ')}。Dockerで起動してください。画面の閲覧は可能です。`:`解析環境は準備済みです · 結果保存先の空き容量 ${health.disk_free_gb} GB · 解析は1件ずつ実行します`;
    $('health').classList.toggle('warning',!!missing.length);
  } catch(err) {$('health').textContent=`環境確認に失敗しました: ${err.message}`;}
  if (/^#[a-f0-9]{32}$/.test(location.hash)) choose(location.hash.slice(1));
  await refresh();setInterval(refresh,3000);
}
init();
