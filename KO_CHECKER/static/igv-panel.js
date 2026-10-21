'use strict';
(() => {
  const panel=document.createElement('section');
  panel.className='card';
  const heading=document.createElement('h2');heading.textContent='IGV：WT・KOの実データ';
  const note=document.createElement('p');note.className='hint';
  note.textContent='カバレッジ、RNA-seqリード、エクソン・イントロン構造を比較します。拡大すると塩基や変異を確認できます。カバレッジは読み取り深度で、ライブラリ量を補正した発現比ではありません。';
  const button=document.createElement('button');button.type='button';button.textContent='IGVを表示';
  const status=document.createElement('p');status.setAttribute('role','status');
  const viewer=document.createElement('div');viewer.id='igv-browser';
  panel.append(heading,note,button,status,viewer);
  document.getElementById('structure-title').closest('section').before(panel);
  let jobId=null, browser=null, generation=0, library=null;
  function loadLibrary(){
    if(globalThis.igv)return Promise.resolve(globalThis.igv);
    if(!library)library=new Promise((resolve,reject)=>{
      const script=document.createElement('script');script.src='/vendor/igv-3.0.0.min.js';
      script.onload=()=>resolve(globalThis.igv);script.onerror=()=>{library=null;script.remove();reject(new Error('IGVライブラリを読み込めませんでした。'));};
      document.head.append(script);
    });
    return library;
  }
  globalThis.prepareIGV=id=>{
    if(jobId===id)return;
    jobId=id;generation++;
    if(browser){globalThis.igv.removeBrowser(browser);browser=null;}
    viewer.replaceChildren();button.hidden=false;button.disabled=false;status.textContent='';
  };
  button.addEventListener('click',async()=>{
    const current=++generation,id=jobId;button.disabled=true;status.textContent='IGVと解析データを読み込み中…';
    let created=null;
    try{
      const response=await fetch(`/api/jobs/${id}/igv/config`);
      const config=await response.json();if(!response.ok)throw new Error(config.error||'表示データがありません。');
      const igv=await loadLibrary();if(current!==generation)return;
      created=await igv.createBrowser(viewer,config);
      if(current!==generation){igv.removeBrowser(created);return;}
      browser=created;button.hidden=true;status.textContent='遺伝子注釈と検体別の読み取りを表示しています。ドラッグで移動、＋／−で拡大・縮小できます。';
    }catch(error){if(current===generation){viewer.replaceChildren();status.textContent=`表示できません：${error.message}`;}}
    finally{if(current===generation)button.disabled=false;}
  });
})();
