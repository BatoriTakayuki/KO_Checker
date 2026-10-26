'use strict';
(() => {
  let entries=[], selection=new Set(), checks=new Map(), token=null;
  let busy=false, blocked=false, limit=300;
  function bytes(value) {
    const amount=Number(value) || 0, units=['B','KB','MB','GB','TB'];
    const scale=amount>0?Math.min(4,Math.floor(Math.log(amount)/Math.log(1024))):0;
    return `${number(amount/1024**scale,scale?2:0)} ${units[scale]}`;
  }
  function buttons() {
    for (const id of ['storage-select-reads','storage-select-bam','storage-select-reference','storage-clear']) $(id).disabled=busy || blocked;
    $('storage-refresh').disabled=busy;
    $('storage-preview').disabled=busy || blocked || !selection.size || selection.size>limit;
    $('storage-delete').disabled=busy || blocked || !token;
    $('storage-cancel').disabled=busy;
    for (const row of entries) {
      const checkbox=checks.get(row.id);
      if (checkbox) checkbox.disabled=busy || blocked || !row.deletable;
    }
  }
  function changed() {
    token=null;$('storage-confirmation').hidden=true;
    const size=entries.filter(e=>selection.has(e.id)).reduce((sum,e)=>sum+e.size_bytes,0);
    $('storage-selected').textContent=selection.size?`${selection.size}件選択 · ${bytes(size)}`:'選択なし';
    $('storage-error').textContent=selection.size>limit?`一度に選択できるのは${limit}件までです。`:'';
    buttons();
  }
  function table(data) {
    entries=data.entries;selection.clear();checks.clear();blocked=data.busy;limit=data.max_selection;
    const target=$('storage-table');target.replaceChildren();
    const head=el('tr');['選択','種類','対象','ファイル容量','削除後の再利用'].forEach(label=>head.append(el('th',label)));
    const thead=el('thead');thead.append(head);target.append(thead);
    const body=el('tbody');
    for (const entry of [...entries].sort((a,b)=>b.size_bytes-a.size_bytes)) {
      const row=el('tr'),choice=el('td'),checkbox=el('input');checkbox.type='checkbox';
      checkbox.setAttribute('aria-label',`${entry.kind_label} ${entry.label}を選択`);
      checkbox.addEventListener('change',()=>{if(checkbox.checked)selection.add(entry.id);else selection.delete(entry.id);changed();});
      checks.set(entry.id,checkbox);choice.append(checkbox);row.append(choice,el('td',entry.kind_label));
      const name=el('td',entry.label);name.append(el('small',entry.path));
      row.append(name,el('td',bytes(entry.size_bytes)),el('td',entry.deletable?entry.description:entry.error));body.append(row);
    }
    if (!entries.length) {const row=el('tr'),cell=el('td','管理対象の保存データはありません。');cell.colSpan=5;row.append(cell);body.append(row);}
    target.append(body);changed();
    $('storage-warning').textContent=data.warnings.join(' ');
    $('storage-status').replaceChildren('管理対象の合計： ', el('strong',bytes(data.total_bytes)),
      '\n保存先の空き容量 ：', el('strong',bytes(data.free_bytes)));
    if (data.busy) $('storage-status').append('\n実行中・待機中の解析があるため削除できません。完了するか、すべて中断してから容量を更新してください。');
  }
  async function load() {
    busy=true;buttons();$('storage-status').textContent='保存データの容量を確認中…';$('storage-error').textContent='';
    try {table(await api('/api/storage'));}
    catch(err) {$('storage-status').textContent='容量を確認できませんでした。';$('storage-error').textContent=err.message;}
    finally {busy=false;buttons();}
  }
  function select(kinds) {
    if (busy || blocked) return;
    for (const entry of entries) {
      if (entry.deletable && kinds.includes(entry.kind)) {selection.add(entry.id);checks.get(entry.id).checked=true;}
    }
    changed();
  }
  $('storage-panel').addEventListener('toggle',()=>{if($('storage-panel').open && !busy)return load();});
  $('storage-refresh').addEventListener('click',load);
  $('storage-select-reads').addEventListener('click',()=>select(['sra','fastq','trimmed']));
  $('storage-select-bam').addEventListener('click',()=>select(['bam','samples']));
  $('storage-select-reference').addEventListener('click',()=>select(['reference','index']));
  $('storage-clear').addEventListener('click',()=>{selection.clear();for(const checkbox of checks.values())checkbox.checked=false;changed();});
  $('storage-cancel').addEventListener('click',()=>{token=null;$('storage-confirmation').hidden=true;buttons();});
  $('storage-preview').addEventListener('click',async()=>{
    if (busy || blocked || !selection.size) return;
    busy=true;token=null;buttons();$('storage-error').textContent='';$('storage-confirmation').hidden=true;
    try {
      const preview=await api('/api/storage/preview',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({ids:[...selection]})});
      token=preview.token;$('storage-confirm-size').textContent=`削除対象: ${preview.entries.length}件 · ${bytes(preview.total_bytes)}（確認の有効期間は5分です）`;
      $('storage-confirm-list').replaceChildren(...preview.entries.map(e=>el('li',`${e.kind_label} / ${e.label} · ${bytes(e.size_bytes)} — ${e.description}`)));
      $('storage-confirmation').hidden=false;
    } catch(err) {$('storage-error').textContent=err.message;}
    finally {busy=false;buttons();}
  });
  $('storage-delete').addEventListener('click',async()=>{
    if (!token || busy || blocked) return;
    busy=true;buttons();$('storage-error').textContent='';
    const confirmedToken=token;token=null;
    try {
      const result=await api('/api/storage/delete',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({token:confirmedToken,confirmed:true})});
      await load();
      $('storage-status').textContent=`${result.deleted.length}件を削除しました。削除したファイル容量: ${bytes(result.deleted_bytes)}。結論・CSV・ログ・履歴は保持しています。`;
      if (result.failed.length) $('storage-error').textContent=result.failed.map(e=>`${e.kind_label} / ${e.label}: ${e.error}`).join('\n');
    } catch(err) {$('storage-error').textContent=err.message;$('storage-confirmation').hidden=true;}
    finally {busy=false;buttons();}
  });
})();
