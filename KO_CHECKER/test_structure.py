import unittest

import rnako
from structure import reference_model, loss_model, variant_model, coordinate_template
from test_core import FakeReference


class StructureTests(unittest.TestCase):
    def setUp(self):
        self.sequence='ATGAAACCC'+'N'*11+'CCCGGGAAA'+'N'*11+'TTTTAACCC'
        self.fa=FakeReference(self.sequence)
        self.fa.references=('chr1',)
        self.tx=rnako.Transcript('T','G','Test','chr1','+','protein_coding',
            [(0,9),(20,29),(40,49)],[(0,9),(20,29),(40,43)])

    def test_reference_codons_and_protein_length(self):
        model=reference_model(self.tx,self.fa)
        self.assertEqual(model['start']['bases'],'ATG')
        self.assertEqual(model['start']['genomic_positions_1based'],[1,2,3])
        self.assertEqual(model['stop']['bases'],'TAA')
        self.assertEqual(model['stop']['genomic_positions_1based'],[44,45,46])
        self.assertEqual(model['protein_length'],7)

    def test_minus_strand_codons_are_mapped_in_transcription_order(self):
        n=len(self.sequence)
        tx=rnako.Transcript('T','G','Test','chr1','-','protein_coding',
            [(n-b,n-a) for a,b in self.tx.exons],[(n-b,n-a) for a,b in self.tx.cds])
        fa=FakeReference(rnako.rc(self.sequence));fa.references=('chr1',)
        model=reference_model(tx,fa)
        self.assertEqual(model['start']['bases'],'ATG')
        self.assertEqual(model['start']['genomic_positions_1based'],[49,48,47])
        self.assertEqual(model['stop']['genomic_positions_1based'],[6,5,4])

    def test_codon_spanning_exons_keeps_all_three_genomic_positions(self):
        tx=rnako.Transcript('T','G','Test','chr1','+','protein_coding',
            [(0,1),(20,28)],[(0,1),(20,25)])
        seq='A'+'N'*19+'TGAAATAA'
        fa=FakeReference(seq);fa.references=('chr1',)
        model=reference_model(tx,fa)
        self.assertEqual(model['start']['bases'],'ATG')
        self.assertEqual(model['start']['genomic_positions_1based'],[1,21,22])

    def test_loss_model_stop_is_projected_back_to_surviving_exon(self):
        row=dict(event_id='loss',lost_exons='2',stop_codon_cdna_1based='13',
            evidence='coverage_only',effect='region_loss_inframe',predicted_protein='MKPF',
            protein_prediction_status='conditional_exon_join_model',ko_support='yes',assumptions='conditional')
        model=loss_model(self.tx,row,self.fa)
        self.assertEqual(model['stop']['bases'],'TAA')
        self.assertEqual(model['stop']['genomic_positions_1based'],[44,45,46])
        self.assertEqual(model['exon_numbers'],[1,3])
        self.assertEqual(model['protein_length'],4)
        self.assertEqual(model['coding_segments'],[{'start_1based':1,'end_1based':9},{'start_1based':41,'end_1based':43}])

    def test_start_or_terminal_loss_does_not_invent_codons(self):
        row=dict(event_id='loss',lost_exons='1',stop_codon_cdna_1based='',
            evidence='coverage_only',effect='start_region_lost',predicted_protein='',
            protein_prediction_status='start_lost',ko_support='yes',assumptions='conditional')
        model=loss_model(self.tx,row,self.fa)
        self.assertIsNone(model['start']);self.assertIsNone(model['stop'])
        row.update(lost_exons='3',effect='terminal_coding_region_loss')
        self.assertIsNone(loss_model(self.tx,row,self.fa)['stop'])

    def test_stop_gained_variant_has_new_stop_and_start_change_is_untranslated(self):
        event=dict(event_id='snv',pos_1based='4',ref='A',alt='T')
        model=variant_model(self.tx,event,self.fa)
        self.assertEqual(model['stop']['bases'],'TAA')
        self.assertEqual(model['stop']['genomic_positions_1based'],[4,5,6])
        self.assertEqual(model['protein_length'],1)
        self.assertEqual(model['coding_segments'],[{'start_1based':1,'end_1based':3}])
        event.update(pos_1based='1',ref='A',alt='C')
        model=variant_model(self.tx,event,self.fa)
        self.assertEqual(model['start']['bases'],'CTG')
        self.assertEqual(model['start']['status'],'start_changed')
        self.assertIsNone(model['stop']);self.assertIsNone(model['protein_length'])

    def test_insertion_keeps_coordinate_mapping_length_consistent(self):
        event=dict(event_id='ins',pos_1based='4',ref='A',alt='AT')
        model=variant_model(self.tx,event,self.fa)
        self.assertEqual(model['transcript_length'],28)
        self.assertEqual(model['start']['genomic_positions_1based'],[1,2,3])
        self.assertEqual(model['effect'],'frameshift')


if __name__=='__main__':unittest.main()
