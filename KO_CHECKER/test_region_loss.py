import csv
import json
from pathlib import Path
import tempfile
import unittest

import rnako
from region_loss import assess_regions, model_loss, write_region_loss, apply_region_classes, REFERENCE_REASON
from reporting import read_conclusion


class Reference:
    def __init__(self, seq): self.seq = seq
    def fetch(self, chrom, a, b): return self.seq[a:b]


class RegionTests(unittest.TestCase):
    def test_reference_region_retention_is_local_and_label_independent(self):
        samples = [dict(srx_id=sid, metadata_label='KO', predicted_class='undetermined',
                        unique_primary_library_fragments=100) for sid in ('retained', 'low', 'missing', 'deleted')]
        exons = [dict(srx_id=sid, chrom='chr1', exon_start_1based=21, exon_end_1based=26,
                      mean_depth=depth, sufficient_depth_fraction=fraction)
                 for sid, depth, fraction in [('retained', 20, 1), ('low', 20, .7)]]
        models = [dict(srx_id='deleted', ko_support='yes', lost_regions='chr1:21-26')]
        apply_region_classes(samples, exons, models)
        self.assertEqual([s['predicted_class'] for s in samples],
                         ['WT_like', 'undetermined', 'undetermined', 'KO_like'])
        self.assertEqual(samples[0]['classification_reason'], REFERENCE_REASON)
        apply_region_classes(samples, exons, [])
        self.assertTrue(all(s['predicted_class'] == 'undetermined' for s in samples))

    def test_reference_retention_preserves_other_evidence_and_requires_all_regions(self):
        sample = dict(srx_id='s', predicted_class='KO_like', classification_reason='sequence event',
                      unique_primary_library_fragments=100)
        models = [dict(srx_id='other', ko_support='yes', lost_regions='chr1:21-26;chr1:41-49')]
        exon = dict(srx_id='s', chrom='chr1', exon_start_1based=21, exon_end_1based=26,
                    mean_depth=100, sufficient_depth_fraction=1)
        apply_region_classes([sample], [exon], models)
        self.assertEqual(sample['classification_reason'], 'sequence event')
        sample['predicted_class'] = 'undetermined'
        apply_region_classes([sample], [exon], models)
        self.assertEqual(sample['predicted_class'], 'undetermined')

    def setUp(self):
        self.tx = rnako.Transcript('T', 'G', 'Test', 'chr1', '+', 'protein_coding',
            [(0, 9), (20, 26), (40, 49)], [(0, 9), (20, 26), (40, 43)])
        self.ref = Reference('ATGAAACCC'+'N'*11+'CCCGGG'+'N'*14+'TTTTAACCC')

    def test_internal_inframe_model_has_correct_peptide(self):
        result = model_loss(self.tx, [(20, 26)], self.ref, 'chr1')
        self.assertEqual(result['reference_protein'], 'MKPPGF')
        self.assertEqual(result['predicted_protein'], 'MKPF')
        self.assertEqual(result['effect'], 'region_loss_inframe')
        self.assertEqual(result['deleted_coding_nt'], 6)

    def test_minus_strand_uses_transcription_order(self):
        seq = rnako.rc(self.ref.seq)
        n = len(seq)
        tx = rnako.Transcript('R', 'G', 'Test', 'chr1', '-', 'protein_coding',
            [(n-b,n-a) for a,b in self.tx.exons], [(n-b,n-a) for a,b in self.tx.cds])
        result = model_loss(tx, [(n-26,n-20)], Reference(seq), 'chr1')
        self.assertEqual(result['predicted_protein'], 'MKPF')
        self.assertEqual(result['reference_protein'], 'MKPPGF')

    def test_start_loss_does_not_invent_alternative_initiation(self):
        result = model_loss(self.tx, [(0, 9)], self.ref, 'chr1')
        self.assertEqual(result['effect'], 'start_region_lost')
        self.assertEqual(result['predicted_protein'], '')

    def test_terminal_loss_is_only_a_known_fragment(self):
        result = model_loss(self.tx, [(40, 49)], self.ref, 'chr1')
        self.assertEqual(result['predicted_protein'], 'MKPPG')
        self.assertIn('unknown_3prime_end', result['protein_prediction_status'])
        self.assertEqual(result['nmd_prediction'], 'unknown')
        self.assertEqual(result['stop_codon_cdna_1based'], '')

    def test_joint_deletion_and_frameshift(self):
        tx = rnako.Transcript('T', 'G', 'Test', 'chr1', '+', 'protein_coding',
            [(0,9),(20,24),(40,49)], [(0,9),(20,24),(40,43)])
        result = model_loss(tx, [(20,24)], self.ref, 'chr1')
        self.assertEqual(result['effect'], 'region_loss_frameshift')
        result = model_loss(self.tx, [(20,26),(40,49)], self.ref, 'chr1')
        self.assertEqual(result['predicted_protein'], 'MKP')
        self.assertEqual(result['deleted_coding_nt'], 9)

    def test_utr_only_loss_keeps_protein(self):
        tx = rnako.Transcript('T','G','Test','chr1','+','protein_coding',
            [(0,9),(20,26),(40,49)], [(40,46)])
        ref = Reference('CCCGGGCCC'+'N'*11+'CCCGGG'+'N'*14+'ATGAAATAA')
        result = model_loss(tx, [(20,26)], ref, 'chr1')
        self.assertEqual(result['effect'], 'noncoding_region_loss')
        self.assertEqual(result['reference_protein'], result['predicted_protein'])

    def data(self):
        samples = [dict(srx_id=sid,metadata_label=label,unique_primary_library_fragments=lib)
                   for sid,label,lib in [('W','WT',1000),('K','KO',2000)]]
        rows=[]
        for sid,label,depths in [('W','WT',[100,100,100]),('K','KO',[200,4,200])]:
            for (a,b),depth in zip(self.tx.exons,depths):
                rows.append(dict(srx_id=sid,metadata_label=label,exon_start_1based=a+1,exon_end_1based=b,
                    mean_depth=depth,sufficient_depth_fraction=1 if depth>=10 else 0,
                    covered_fraction=1,other_exons_sufficient_fraction=0))
        return samples,rows

    def test_residual_reads_and_bad_global_coverage_do_not_hide_loss(self):
        samples,rows=self.data()
        losses,anchors=assess_regions(samples,rows)
        self.assertAlmostEqual(losses['K'][(20,26)][0], .02)
        self.assertEqual(len(anchors['K']),2)

    def test_global_rna_loss_and_low_library_cannot_call_region_loss(self):
        samples,rows=self.data()
        for row in rows:
            if row['srx_id']=='K': row.update(mean_depth=1,sufficient_depth_fraction=0)
        self.assertEqual(assess_regions(samples,rows)[0],{})
        samples,rows=self.data();samples[1]['unique_primary_library_fragments']=0
        self.assertEqual(assess_regions(samples,rows)[0],{})

    def test_low_wt_depth_is_not_evidence_for_loss(self):
        samples,rows=self.data()
        rows[1].update(mean_depth=4,sufficient_depth_fraction=0)
        self.assertNotIn((20,26),assess_regions(samples,rows)[0]['K'])

    @unittest.skipIf(rnako.pysam is None,'pysam required')
    def test_major_coding_loss_updates_verdict_and_exports_conditional_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            seq='ATGAAACCC'+'N'*11+'CCC'*20+'N'*20+'TTTTAACCC'
            fasta=root/'reference.fa';fasta.write_text('>chr1\n'+seq+'\n')
            rnako.pysam.faidx(str(fasta))
            gtf=root/'reference.gtf'
            attrs='gene_id "G"; gene_name "Test"; transcript_id "T"; gene_type "protein_coding";'
            with gtf.open('w') as f:
                for a,b in [(0,9),(20,80),(100,109)]:
                    f.write(f'chr1\ttest\texon\t{a+1}\t{b}\t.\t+\t.\t{attrs}\n')
                    f.write(f'chr1\ttest\tCDS\t{a+1}\t{min(b,103)}\t.\t+\t0\t{attrs}\n')
                # A longer, unresolved CDS (starts TGA) must not mask the
                # canonical-start transcript, as happened with Eed annotations.
                longer=attrs.replace('transcript_id "T"','transcript_id "L"')
                f.write(f'chr1\ttest\texon\t1\t109\t.\t+\t.\t{longer}\n')
                f.write(f'chr1\ttest\tCDS\t2\t109\t.\t+\t0\t{longer}\n')
            samples,rows=self.data()
            for s in samples:
                s.update(predicted_class='undetermined',classification_reason='',gene='Test',
                         gene_unique_fragment_cpm=100,coverage_warning='')
            for r in rows:
                if int(r['exon_start_1based'])==21:r['exon_end_1based']=80
                if int(r['exon_start_1based'])==41:r.update(exon_start_1based=101,exon_end_1based=109)
            for name,data in [('samples',samples),('exon_coverage',rows)]:
                with (root/(name+'.csv')).open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
            (root/'junctions.csv').write_text('srx_id,support_fragments,intron_start_1based,intron_end_1based,comparison_status\n')
            (root/'run.json').write_text(json.dumps({'thresholds':{'min_event_depth':10}}))
            models=write_region_loss(root,fasta,gtf,'Test')
            self.assertEqual(len(models),1)
            self.assertEqual(models[0]['ko_support'],'yes')
            self.assertEqual(models[0]['evidence'],'coverage_only')
            self.assertEqual(models[0]['predicted_protein'],'MKPF')
            self.assertIn('MKPF',(root/'region_proteins.fasta').read_text())
            report=read_conclusion(root)
            self.assertEqual(report['sequence']['status'],'ko_supported')
            self.assertIn('RNA領域欠損',report['sequence']['text'])
            self.assertEqual(report['region_loss']['ko_samples'],['K'])
            self.assertEqual(report['nmd']['status'],'unavailable')
            # A real, enriched junction changes the evidence label, not the
            # assumption that this is a single transcript model rather than DNA.
            (root/'junctions.csv').write_text('srx_id,support_fragments,intron_start_1based,intron_end_1based,comparison_status\nK,12,10,100,enriched_vs_comparators\n')
            self.assertEqual(write_region_loss(root,fasta,gtf,'Test')[0]['evidence'],'junction_supported')
            # A shared/rare junction must not be advertised as KO-specific.
            (root/'junctions.csv').write_text('srx_id,support_fragments,intron_start_1based,intron_end_1based,comparison_status\nK,12,10,100,shared_or_not_enriched\n')
            self.assertEqual(write_region_loss(root,fasta,gtf,'Test')[0]['evidence'],'coverage_only')
            # Repeating this post-processing must not accumulate duplicate models.
            self.assertEqual(len(write_region_loss(root,fasta,gtf,'Test')),1)
            # Removing an old region call when WT is no longer assessable must
            # not leave a stale KO verdict in samples.csv.
            with (root/'exon_coverage.csv').open() as f: reread=list(csv.DictReader(f))
            for row in reread:
                if row['metadata_label']=='WT': row['sufficient_depth_fraction']='0'
            with (root/'exon_coverage.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(reread[0]));w.writeheader();w.writerows(reread)
            self.assertEqual(write_region_loss(root,fasta,gtf,'Test'),[])
            self.assertEqual(read_conclusion(root)['sequence']['counts']['KO_like'],0)


if __name__=='__main__': unittest.main()
