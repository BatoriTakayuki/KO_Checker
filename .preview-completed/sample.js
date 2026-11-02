// Fictional data for an isolated screenshot. No network or analysis requests.
const samples = ['DEMO_WT1','DEMO_WT2','DEMO_KO1','DEMO_KO2'].map((id,i)=>({
  srx_id:id,metadata_label:i<2?'WT':'KO',predicted_class:i<2?'WT_like':'KO_like',
  gene_unique_fragment_cpm:[42.8,45.2,12.4,10.8][i],expression_relative_to_wt:[.973,1.027,.282,.245][i],
  coverage_warning:'',classification_reason:i<2?'候補領域の読み取りを保持（架空例）':'RNA領域欠損を支持（架空例）'
}));
const exons = [1,2,3,4].map((n)=>({number:n,start_1based:100000+(n-1)*1200,end_1based:100299+(n-1)*1200,cds:[{start_1based:100000+(n-1)*1200,end_1based:100299+(n-1)*1200}]}));
const coverage = samples.flatMap((s,i)=>exons.map((e,j)=>({
  srx_id:s.srx_id,metadata_label:s.metadata_label,chrom:'chr1',exon_start_1based:e.start_1based,exon_end_1based:e.end_1based,
  mean_depth:i<2?48+j*3:j===1?1.2:16+j,min_depth:i<2?25:j===1?0:8,overlapping_fragments:i<2?95:24,
  covered_fraction:i>=2&&j===1?.12:1,sufficient_depth_fraction:i>=2&&j===1?.02:.96,
  coverage_status:'covered',comparison_status:i<2?'reference_sample':j===1?'suspected_local_depletion':'no_specific_depletion',
  region_loss_status:i>=2&&j===1?'suspected_region_loss':'',region_normalized_depth_ratio:i<2?1:j===1?.04:.95,
  relative_depth_vs_wt:i<2?1:j===1?.04:.95,wt_reference_samples:2,skipping_junction_fragments:i<2?0:18,
  longest_zero_depth_region:i>=2&&j===1?'chr1:101230-101410':'',transcript_exons:'DEMO_TX1:exon'+(j+1)
})));
const marker=(pos,bases)=>({bases,genomic_positions_1based:[pos,pos+1,pos+2],cdna_1based:1,status:'架空例'});
const reference={protein_length:399,start:marker(100000,'ATG'),stop:marker(103897,'TAA')};
const lossModel={kind:'region_loss',removed_exon_numbers:[2],protein_length:299,evidence:'junction_supported',start:reference.start,stop:reference.stop,note:'架空の表示例：エクソン2の領域欠損を仮定したモデルです。'};
const structure={gene:'DemoGene',samples:samples.map(s=>({id:s.srx_id,group:s.metadata_label})),coverage,variants:[],
  transcripts:[{id:'DEMO_TX1',chrom:'chr1',strand:'+',exons,reference,variant_models:[],loss_models:{DEMO_KO1:lossModel,DEMO_KO2:lossModel}}]};
const regionRows=samples.slice(2).map(s=>({srx_id:s.srx_id,transcript_id:'DEMO_TX1',lost_exons:'2',lost_regions:'chr1:101200-101499',maximum_normalized_depth_ratio:.04,evidence:'junction_supported',junction_support:18,ko_support:'yes',effect:'region_loss_inframe',deleted_coding_nt:300,coding_loss_fraction:.25,protein_prediction_status:'conditional_exon_join_model',nmd_prediction:'表示例：早期終止なし',assumptions:'架空データによる条件付きモデル'}));
const report={sequence:{text:'KO群2検体でRNA領域欠損を支持する所見があります（表示例）。',note:'このプレビューは架空のデータです。実際の解析結果ではありません。'},
 expression:{text:'KO群の平均CPMはWT群の26.4%（約73.6%減少）です。',test:{text:'有意差の検定は未実施です（サンプル表示）。'}},
 exons:{text:'KO群でエクソン2の局所的な読み取り低下を検出した場合の表示例です。',strongest:coverage.filter(r=>r.region_loss_status)},
 nmd:{text:'この仮定モデルでは早期終止はなく、NMDを支持する予測はありません（表示例）。'},
 overall:'表示例：RNA上でKO様の所見を認めます。DNA上の欠失やタンパク質消失を確定するものではありません。',notes:['数値・検体ID・遺伝子構造はすべて架空です。']};
const demoData={tables:{samples,exon_coverage:coverage,region_deletions:regionRows,observations:[],events:[],junctions:[],consequences:[]},conclusion:report,structure,counts:{samples:4,exon_coverage:16,region_deletions:2,events:0,junctions:0},region_analysis_available:true};
$('empty').hidden=true;$('job').hidden=false;
$('health').textContent='完了画面の見本画像 ｜ 数値・判定・検体IDは架空です。実際の解析とは関係ありません。';
$('job-id').textContent='PREVIEW · 架空のサンプル';$('job-title').textContent='DemoGene / mm10';
$('status').textContent='完了';$('status').className='badge completed';
$('job-message').textContent='解析が完了しました。結果を確認してください。';$('progress').value=100;
$('stage').textContent='完了';$('percent').textContent='100%（工程の目安）';
$('elapsed').textContent='今回の実行時間: 3時間24分18秒（架空例）';$('run-count').textContent='run完了: 4 / 4';
$('pause').hidden=true;$('resume').hidden=true;$('retry').hidden=false;
$('gene').value='DemoGene';$('wt').value='DEMO_WT1\nDEMO_WT2';$('ko').value='DEMO_KO1\nDEMO_KO2';
$('history').replaceChildren(el('p','DemoGene · mm10 / 完了（表示例）'));
results(demoData,{id:'preview'});
document.querySelectorAll('a').forEach(node=>node.removeAttribute('href'));
document.title='KO_Checker 完了画面プレビュー（架空データ）';
