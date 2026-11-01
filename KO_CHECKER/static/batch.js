'use strict';
(() => {
  const batchStatus={running:'実行中',paused:'一時停止中',completed:'完了',cancelled:'取り消し済み'};
  const itemStatus={pending:'未実行',queued:'待機中',running:'解析中',paused:'中断中',completed:'完了',failed:'失敗',cancelled:'取り消し'};
  let csvText='', csvName='', rows=[], selection=new Set(), checks=new Map(), selected=null, timer=null;
  const post=(path,body)=>api(path,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
  const gb=value=>value===null||value===undefined?'—':`${number(value,value<1?3:1)} GB`;
  function genes(row) {
    const mode=$('batch-multi-gene').value;
    return row.genes.length<2||mode==='split'?row.genes:mode==='first'?row.genes.slice(0,1):[];
  }
  function changed() {
    const chosen=rows.filter(r=>selection.has(r.index));
    const analyses=chosen.reduce((sum,r)=>sum+genes(r).length,0);
    const skipped=chosen.filter(r=>!genes(r).length).length;
    const fastq=chosen.filter(r=>genes(r).length).reduce((sum,r)=>sum+(r.fastq_gb||0),0);
    $('batch-selected').textContent=chosen.length?`${chosen.length}実験を選択 · ${analyses}解析${skipped?`（複数遺伝子の${skipped}実験はスキップ）`:''} · FASTQ合計 ${gb(fastq)}`:'選択なし';
    $('batch-start').disabled=!analyses;
  }
  function select(predicate) {
    selection=new Set(rows.filter(r=>!r.error&&predicate(r)).map(r=>r.index));
    for (const [index,box] of checks) box.checked=selection.has(index);
    changed();
  }
  function preview(data) {
    rows=data.rows;checks.clear();
    const target=$('batch-preview-table');target.replaceChildren();
    const head=el('tr');['','実験ID','対象遺伝子','WT群のSRX','KO群のSRX','FASTQ合計','CSVの状態','説明','入力エラー'].forEach(t=>head.append(el('th',t)));
    const thead=el('thead');thead.append(head);
    const body=el('tbody');
    for (const row of rows) {
      const tr=el('tr'),box=document.createElement('input');box.type='checkbox';box.disabled=!!row.error;
      box.setAttribute('aria-label',`${row.experiment_id}を選択`);
      box.addEventListener('change',()=>{box.checked?selection.add(row.index):selection.delete(row.index);changed();});
      checks.set(row.index,box);
      const cell=el('td');cell.append(box);tr.append(cell);
      for (const value of [row.experiment_id,row.genes.join(' / '),row.wt,row.ko,gb(row.fastq_gb),row.status,row.description,row.error]) tr.append(el('td',value||'—'));
      if (row.error) tr.className='batch-invalid';
      body.append(tr);
    }
    target.append(thead,body);
    const invalid=rows.filter(r=>r.error).length;
    $('batch-preview-status').textContent=`${csvName}: ${rows.length}行を読み込みました${invalid?`（入力エラー ${invalid}行は選択できません）`:''}。保存先の空き容量: ${number(data.free_gb,1)} GB。\nCSVの状態に「未公開」を含む行は初期状態で選択していません。`;
    $('batch-preview').hidden=false;
    select(r=>!/未公開/.test(r.status));
  }
  $('batch-file').addEventListener('change',async event=>{
    const file=event.target.files[0];$('batch-error').textContent='';$('batch-preview').hidden=true;rows=[];selection.clear();changed();
    if (!file) return;
    try {
      csvText=await file.text();csvName=file.name;
      preview(await post('/api/batches/preview',{csv:csvText}));
    } catch(err) {$('batch-error').textContent=err.message;}
  });
  $('batch-multi-gene').addEventListener('change',changed);
  $('batch-select-all').addEventListener('click',()=>select(()=>true));
  $('batch-clear').addEventListener('click',()=>select(()=>false));
  $('batch-select-size').addEventListener('click',()=>{
    const limit=Number($('batch-size-limit').value);
    if (!(limit>0)) {$('batch-error').textContent='FASTQ合計の上限をGBで入力してください。';return;}
    $('batch-error').textContent='';select(r=>r.fastq_gb!==null&&r.fastq_gb<=limit&&!/未公開/.test(r.status));
  });
  $('batch-start').addEventListener('click',async()=>{
    $('batch-start').disabled=true;$('batch-error').textContent='';
    try {
      const record=await post('/api/batches',{csv:csvText,name:csvName,selected:[...selection],assembly:$('batch-assembly').value,
        threads:$('batch-threads').value,min_free_gb:$('batch-min-free').value,multi_gene:$('batch-multi-gene').value,
        confirm_replicates:$('batch-confirm').checked,ram_work:$('batch-ram').checked,parallel:$('batch-parallel').value});
      $('batch-file').value='';$('batch-preview').hidden=true;rows=[];selection.clear();changed();
      selected=record.id;await refresh();
    } catch(err) {$('batch-error').textContent=err.message;changed();}
  });
  function counts(record) {
    return ['completed','failed','running','queued','paused','pending','cancelled'].filter(k=>record.counts[k])
      .map(k=>`${itemStatus[k]} ${record.counts[k]}`).join(' / ');
  }
  async function control(operation) {
    if (!selected) return;
    $('batch-error').textContent='';
    try {await post(`/api/batches/${selected}/${operation}`,{});await refresh();}
    catch(err) {$('batch-error').textContent=err.message;}
  }
  for (const operation of ['pause','resume','cancel','retry']) $(`batch-${operation}`).addEventListener('click',()=>control(operation));
  function detail(record) {
    $('batch-detail').hidden=false;
    $('batch-title').textContent=`${record.name} · ${record.options.assembly} · ${record.total}解析${record.options.ram_work?' · RAM上で処理':''} · 同時${record.options.parallel||1}件`;
    $('batch-status').textContent=batchStatus[record.status]||record.status;$('batch-status').className=`badge ${record.status}`;
    $('batch-message').textContent=record.message;
    const finished=(record.counts.completed||0)+(record.counts.failed||0);
    $('batch-progress').value=record.total?finished/record.total*100:0;
    $('batch-counts').textContent=counts(record);
    $('batch-pause').hidden=record.status!=='running';
    $('batch-resume').hidden=record.status!=='paused';
    $('batch-cancel').hidden=!['running','paused'].includes(record.status);
    $('batch-retry').hidden=!record.counts.failed;
    $('batch-summary').href=`/api/batches/${record.id}/summary.csv`;
    const target=$('batch-items');target.replaceChildren();
    const head=el('tr');['実験ID','対象遺伝子','状態','工程','メッセージ','結果'].forEach(t=>head.append(el('th',t)));
    const thead=el('thead');thead.append(head);
    const body=el('tbody');
    for (const item of record.items) {
      const tr=el('tr');
      tr.append(el('td',item.experiment_id),el('td',item.gene),el('td',itemStatus[item.status]||item.status),
        el('td',item.job_id&&['running','queued'].includes(item.status)?`${stages[item.stage]||item.stage} ${item.progress}%`:'—'),el('td',item.message||'—'));
      const open=el('td');
      if (item.job_id) {
        const button=el('button','開く','secondary');button.type='button';
        button.addEventListener('click',()=>{choose(item.job_id);$('job').scrollIntoView({behavior:'smooth'});});
        open.append(button);
      } else open.textContent='—';
      tr.append(open);body.append(tr);
    }
    target.append(thead,body);
  }
  async function refresh() {
    try {
      const list=await api('/api/batches');
      const target=$('batch-list');target.replaceChildren();
      for (const record of list) {
        const button=el('button',`${record.name} · ${batchStatus[record.status]||record.status} · ${counts(record)}`);
        button.append(el('small',new Date(record.created_at).toLocaleString('ja-JP')));
        button.className=record.id===selected?'active':'';
        button.addEventListener('click',()=>{selected=record.id;refresh();});
        target.append(button);
      }
      if (!list.length) target.append(el('p','まだ一括処理はありません。','muted'));
      if (!selected && list.length) selected=list[0].id;
      if (selected) detail(await api(`/api/batches/${selected}`));
    } catch(err) {$('batch-error').textContent=`通信エラー: ${err.message}`;}
  }
  $('batch-panel').addEventListener('toggle',()=>{
    clearInterval(timer);
    if ($('batch-panel').open) {refresh();timer=setInterval(refresh,5000);}
  });
})();
