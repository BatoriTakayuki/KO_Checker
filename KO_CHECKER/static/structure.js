'use strict';
(() => {
  const ns='http://www.w3.org/2000/svg',get=id=>document.getElementById(id);
  const palette={wt:'#256aa0',ko:'#a65b28',utr:'#dce9f3',loss:'#b32c3b',start:'#117044',stop:'#793c9e',variant:'#bb6500',ink:'#20352f'};
  let data=null,currentSVG=null;
  function svg(tag,attrs={},content) {
    const node=document.createElementNS(ns,tag);
    for(const [key,value] of Object.entries(attrs))node.setAttribute(key,String(value));
    if(content!==undefined)node.textContent=String(content);
    return node;
  }
  function line(root,x1,y1,x2,y2,attrs={}) {root.append(svg('line',{x1,y1,x2,y2,stroke:'#a6b5b0','stroke-width':1.4,...attrs}));}
  function text(root,x,y,value,attrs={}) {root.append(svg('text',{x,y,fill:palette.ink,'font-size':12,...attrs},value));}
  function interactive(node,title,detail,select) {
    node.append(svg('title',{},title));node.setAttribute('tabindex','0');node.setAttribute('role','button');node.setAttribute('aria-label',title);
    node.style.cursor='pointer';
    const show=()=>select(detail);
    node.addEventListener('click',show);
    node.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();show();}});
  }
  function layout(tx,compressed) {
    const intervals=tx.exons, left=95, span=1100;
    const minimum=Math.min(...intervals.map(e=>e.start_1based)),maximum=Math.max(...intervals.map(e=>e.end_1based));
    const total=intervals.reduce((n,e)=>n+e.end_1based-e.start_1based+1,0);
    let cursor=left;
    const boxes=intervals.map(exon=>{
      let x,width;
      if(compressed) {
        x=cursor;width=Math.max(38,(span-(intervals.length-1)*28)*(exon.end_1based-exon.start_1based+1)/total);cursor+=width+28;
      }else{
        const a=tx.strand==='-'?maximum-exon.end_1based:exon.start_1based-minimum;
        x=left+a/Math.max(1,maximum-minimum)*span;
        width=Math.max(.8,(exon.end_1based-exon.start_1based+1)/Math.max(1,maximum-minimum)*span);
      }
      return {...exon,x,width};
    });
    const map=(position)=>{
      const exon=boxes.find(e=>e.start_1based<=position&&position<=e.end_1based);
      if(!exon)return null;
      const fraction=tx.strand==='-'?(exon.end_1based-position)/(exon.end_1based-exon.start_1based+1):(position-exon.start_1based)/(exon.end_1based-exon.start_1based+1);
      return exon.x+fraction*exon.width;
    };
    return {boxes,map,width:Math.max(1280,...boxes.map(e=>e.x+e.width+65)),minimum,maximum};
  }
  function coverage(source,sample,exon) {
    return source.coverage.find(r=>r.srx_id===sample&&Number(r.exon_start_1based)===exon.start_1based&&Number(r.exon_end_1based)===exon.end_1based);
  }
  function draw(source,tx,wt,ko,model,compressed=true,select=()=>{}) {
    const axis=layout(tx,compressed),height=710;
    const root=svg('svg',{xmlns:ns,viewBox:`0 0 ${axis.width} ${height}`,width:axis.width,height,
      role:'img','aria-label':`${source.gene} ${tx.id} WTとKOのエクソン・イントロン構造`,
      'font-family':'Segoe UI, Yu Gothic, Meiryo, sans-serif'});
    root.append(svg('rect',{width:axis.width,height,fill:'#fff'}));
    const defs=svg('defs'),pattern=svg('pattern',{id:'region-loss-hatch',width:6,height:6,patternUnits:'userSpaceOnUse'});
    pattern.append(svg('rect',{width:6,height:6,fill:'#fff1f1'}),svg('path',{d:'M-1 1L5 7M1 -1L7 5',stroke:'#d17880','stroke-width':1}));defs.append(pattern);root.append(defs);
    text(root,24,29,`${source.gene}  ·  ${tx.id}`,{'font-size':19,'font-weight':700});
    text(root,24,53,`${tx.chrom}:${axis.minimum.toLocaleString()}–${axis.maximum.toLocaleString()}  /  ${tx.strand}鎖  /  左5′ → 右3′${tx.strand==='-'?'（ゲノム座標は右へ減少）':''}`);
    text(root,24,76,compressed?'模式図：イントロン短縮。小さいエクソンは視認性のため幅を拡大。':'ゲノム座標に比例した表示。細かい領域は倍率を上げて確認。',{'font-size':11,fill:'#61736c'});
    const legend=[['#dce9f3','非CDS / UTR'],[palette.wt,'CDS / 予測ORF'],['url(#region-loss-hatch)','欠損候補'],[palette.start,'開始'],[palette.stop,'終止'],[palette.variant,'観測変異']];
    legend.forEach(([color,label],i)=>{const x=24+i*160;root.append(svg('rect',{x,y:92,width:15,height:12,fill:color,stroke:'#8a9993'}));text(root,x+23,103,label,{'font-size':11});});
    text(root,24,129,model?`KO：${model.kind==='region_loss'?'領域欠損の仮定モデル':'1変異のみを適用した仮定モデル'}。DNA欠失・実際の転写産物を確定しません。`:'KO：参照の構造・コドン位置に観測値を重ねています。正常な全長RNAの存在を確定しません。',{'font-size':11,fill:'#85512f'});
    if(model?.kind==='region_loss')text(root,24,148,model.evidence==='junction_supported'?'接合の増加を支持するRNA読み取りあり。':'接合未確認：破線の接合と予測終止コドンは仮定に基づきます。',{'font-size':11,fill:palette.loss});
    function track(sample,label,y,color,trackModel) {
      text(root,24,y-13,label,{'font-size':16,'font-weight':700,fill:color});text(root,24,y+8,sample||'検体なし',{'font-size':10});
      const removed=new Set(trackModel?.removed_exon_numbers||[]);
      for(let i=0;i<axis.boxes.length-1;i++){
        const a=axis.boxes[i],b=axis.boxes[i+1],conditional=!!trackModel;
        line(root,a.x+a.width,y,b.x,y,{stroke:color,'stroke-dasharray':conditional?'5 4':'none'});
        if(b.x-a.x-a.width>12)text(root,(a.x+a.width+b.x)/2,y-4,'›',{'text-anchor':'middle',fill:color});
      }
      // Connect surviving exons across lost blocks; keep absent exons visible
      // underneath so WT and KO always share the same reference-coordinate axis.
      const retained=axis.boxes.filter(e=>!removed.has(e.number));
      for(let i=0;i<retained.length-1;i++){
        const a=retained[i],b=retained[i+1];
        if(b.number-a.number>1)root.append(svg('path',{d:`M${a.x+a.width} ${y-8} Q${(a.x+a.width+b.x)/2} ${y-52} ${b.x} ${y-8}`,
          fill:'none',stroke:color,'stroke-width':2,'stroke-dasharray':trackModel.evidence==='junction_supported'?'none':'6 5'}));
      }
      for(const exon of axis.boxes){
        const loss=removed.has(exon.number),row=coverage(source,sample,exon),group=svg('g');
        group.append(svg('rect',{x:exon.x,y:y-7,width:exon.width,height:14,rx:1,fill:loss?'url(#region-loss-hatch)':palette.utr,stroke:loss?palette.loss:color,'stroke-dasharray':loss?'4 3':'none'}));
        const coding=(!loss&&Array.isArray(trackModel?.coding_segments))?trackModel.coding_segments
          .map(c=>({start_1based:Math.max(c.start_1based,exon.start_1based),end_1based:Math.min(c.end_1based,exon.end_1based)}))
          .filter(c=>c.start_1based<=c.end_1based):exon.cds;
        for(const cds of coding){
          const a=axis.map(cds.start_1based),b=axis.map(cds.end_1based);
          group.append(svg('rect',{x:Math.min(a,b),y:y-12,width:Math.max(1,Math.abs(a-b)+exon.width/(exon.end_1based-exon.start_1based+1)),height:24,
            fill:loss?'url(#region-loss-hatch)':color,stroke:loss?palette.loss:color,'stroke-dasharray':loss?'4 3':'none'}));
        }
        if(exon.width>=23)text(group,exon.x+exon.width/2,y+29,`E${exon.number}`,{'text-anchor':'middle','font-size':10,fill:loss?palette.loss:palette.ink});
        if(exon.width>=32)text(group,exon.x+exon.width/2,y+45,row?`${Number(row.mean_depth).toFixed(1)}×`:'未計測',{'text-anchor':'middle','font-size':10,fill:'#61736c'});
        const detail=`${label} ${sample} · E${exon.number} · ${tx.chrom}:${exon.start_1based}-${exon.end_1based} · 平均深度 ${row?.mean_depth??'未計測'} · 読まれた塩基 ${row?`${(Number(row.covered_fraction)*100).toFixed(1)}%`:'未計測'}${loss?' · 欠損候補（RNA上の推定）':''}`;
        interactive(group,detail,detail,select);root.append(group);
      }
      const markers=trackModel||tx.reference;
      for(const [kind,dy,colorMark] of [['start',-45,palette.start],['stop',-75,palette.stop]]){
        const marker=markers[kind];if(!marker)continue;
        const positions=marker.genomic_positions_1based.map(axis.map).filter(x=>x!==null);if(!positions.length)continue;
        const x=positions[0],g=svg('g');
        line(g,x,y-13,x,y+dy+5,{stroke:colorMark,'stroke-dasharray':trackModel?'3 3':'none'});
        g.append(kind==='start'?svg('path',{d:`M${x} ${y+dy-5} l-6 10 h12 Z`,fill:colorMark}):svg('rect',{x:x-5,y:y+dy-5,width:10,height:10,fill:colorMark}));
        const startUnresolved=kind==='start'&&marker.bases!=='ATG';
        text(g,Math.min(axis.width-125,Math.max(80,x+9)),y+dy,`${kind==='start'?'開始':'終止'} ${marker.bases}${startUnresolved?'（開始未確認）':trackModel?'（仮定）':''}`,{'font-size':10,fill:colorMark});
        const detail=`${label} ${kind==='start'?'開始':'終止'}コドン ${marker.bases} · cDNA ${marker.cdna_1based}–${marker.cdna_1based+2} · ${tx.chrom} ゲノム位置 ${marker.genomic_positions_1based.join(', ')} · ${marker.status}。${trackModel?'条件付きモデル。':'参照注釈上の位置で、WT配列の実測確定ではありません。'}`;
        interactive(g,detail,detail,select);root.append(g);
      }
      const variantGroups=new Map();
      for(const variant of source.variants.filter(v=>v.srx_id===sample&&Number(v.alt_reads)>0)){
        if(!variantGroups.has(variant.pos_1based))variantGroups.set(variant.pos_1based,[]);
        variantGroups.get(variant.pos_1based).push(variant);
      }
      for(const variants of variantGroups.values()){
        const variant=variants[0],x=axis.map(Number(variant.pos_1based));if(x===null)continue;
        const g=svg('g'),vy=y+65;
        g.append(svg('path',{d:`M${x} ${vy-6} l6 6 -6 6 -6 -6 Z`,fill:palette.variant,stroke:'#fff'}));
        const short=variants.length>1?`同位置に変異${variants.length}種`:variant.ref.length+variant.alt.length<12?`${variant.ref}→${variant.alt} · ${variant.alt_reads}/${variant.depth_at_anchor}`:`挿入欠失 · ${variant.alt_reads}/${variant.depth_at_anchor}`;
        text(g,x,vy+17,short,{'text-anchor':'middle','font-size':10,fill:palette.variant});
        const detail=variants.map(v=>`${label} ${sample} · ${v.chrom}:${v.pos_1based} ${v.ref}→${v.alt} · 変異支持 ${v.alt_reads}/${v.depth_at_anchor} · ${(Number(v.local_alt_fraction)*100).toFixed(1)}%${v.coverage_status==='low_depth'?'（低深度・参考値）':''}`).join('\n')+'。検出した変異がKOの原因とは限りません。';
        interactive(g,detail,detail,select);root.append(g);
      }
    }
    track(wt,'WT',255,palette.wt,null);track(ko,'KO',450,palette.ko,model);
    text(root,24,562,'箱の下：エクソン番号と観測平均深度。変異・箱・コドンをクリック／Enterで詳細。',{'font-size':11,fill:'#61736c'});
    const maxProtein=Math.max(tx.reference.protein_length||1,model?.protein_length||1),barWidth=560;
    text(root,24,595,'予測タンパク質長',{'font-size':12,'font-weight':600});
    for(const [value,label,y,color] of [[tx.reference.protein_length,'参照',616,palette.wt],[model?.protein_length,'KOモデル',646,palette.ko]]){
      text(root,24,y+9,label,{'font-size':11});
      if(value!==null&&value!==undefined){root.append(svg('rect',{x:110,y,width:value/maxProtein*barWidth,height:12,fill:color}));text(root,125+value/maxProtein*barWidth,y+10,`${value} aa${label==='KOモデル'?'（条件付き）':''}`,{'font-size':11});}
      else text(root,110,y+10,'予測なし／翻訳開始・末端などが未確定',{'font-size':11});
    }
    text(root,axis.width-24,690,'参照注釈＋RNA観測／条件付きモデル。DNA上のKOを確定する図ではありません。',{'text-anchor':'end','font-size':10,fill:'#61736c'});
    return root;
  }
  function options(id,items,value) {
    const node=get(id);node.replaceChildren();
    for(const [key,label] of items){const option=document.createElement('option');option.value=key;option.textContent=label;node.append(option);}
    if(items.some(([key])=>key===value))node.value=value;
  }
  function selectedTranscript(){return data?.transcripts.find(t=>t.id===get('structure-transcript').value);}
  function models() {
    const tx=selectedTranscript();if(!tx)return;
    const sid=get('structure-ko').value,previous=get('structure-model').value;
    const choices=[['reference','参照位置＋KO観測値（モデル適用なし）']];
    if(tx.loss_models[sid])choices.push(['region','領域欠損を仮定したモデル']);
    for(const model of tx.variant_models.filter(m=>m.srx_id===sid))choices.push([model.event_id,`1変異モデル ${model.event_id} · ${model.effect}`]);
    options('structure-model',choices,choices.some(([key])=>key===previous)?previous:tx.loss_models[sid]?'region':'reference');
  }
  function redraw(){
    const tx=selectedTranscript();if(!tx)return;
    const sid=get('structure-ko').value,key=get('structure-model').value;
    const model=key==='region'?tx.loss_models[sid]:tx.variant_models.find(m=>m.srx_id===sid&&m.event_id===key);
    currentSVG=draw(data,tx,get('structure-wt').value,sid,model, get('structure-compress').checked,
      detail=>{get('structure-detail').textContent=detail;});
    currentSVG.style.width=`${Number(get('structure-zoom').value||1)*100}%`;currentSVG.style.height='auto';
    get('structure-diagram').replaceChildren(currentSVG);
    get('structure-detail').textContent=model?.note||'エクソン・開始／終止コドン・変異を選ぶと詳細を表示します。';
  }
  globalThis.renderGeneStructure=source=>{
    data=source;currentSVG=null;get('structure-diagram').replaceChildren();
    get('structure-controls').hidden=!source?.transcripts?.length;
    get('structure-save').disabled=!source?.transcripts?.length;
    get('structure-status').textContent=source?.transcripts?.length?'転写産物と検体を選んで比較できます。接合が未確認のKOモデルは仮定として表示します。':'この結果には構造図データがありません。同じ条件で再実行すると追加できます。';
    if(!source?.transcripts?.length){get('structure-detail').textContent='';return;}
    options('structure-transcript',source.transcripts.map(t=>[t.id,`${t.id} · ${t.reference.protein_length??'未確定'} aa`]),source.transcripts[0].id);
    for(const [id,group] of [['structure-wt','WT'],['structure-ko','KO']]){
      const samples=source.samples.filter(s=>s.group===group);options(id,samples.map(s=>[s.id,s.id]),samples[0]?.id);
    }
    // A new result should default to its region-loss model, never retain a
    // reference-only selection made for an unrelated previous analysis.
    get('structure-model').replaceChildren();get('structure-model').value='';models();redraw();
  };
  for(const id of ['structure-transcript','structure-ko'])get(id).addEventListener('change',()=>{models();redraw();});
  for(const id of ['structure-wt','structure-model','structure-compress','structure-zoom'])get(id).addEventListener('change',redraw);
  get('structure-save').addEventListener('click',()=>{
    if(!currentSVG)return;
    const copy=currentSVG.cloneNode(true);copy.removeAttribute('style');
    const blob=new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml;charset=utf-8'}),url=URL.createObjectURL(blob);
    const link=document.createElement('a');link.href=url;link.download=`${data.gene}_${selectedTranscript().id}_WT_KO.svg`;link.click();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  globalThis.GeneStructureDiagram={draw,layout};
})();
