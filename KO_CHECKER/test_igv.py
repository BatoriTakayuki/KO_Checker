import json
from pathlib import Path
import tempfile
import unittest

from igv_view import byte_range, confined_file, annotations


class IGVTests(unittest.TestCase):
    def test_byte_ranges(self):
        for header, expected in [('bytes=0-9',(0,9)),('bytes=5-',(5,99)),
                                 ('bytes=-10',(90,99)),('bytes=90-200',(90,99))]:
            self.assertEqual(byte_range(header,100),expected)
        for header in ('bytes=100-', 'bytes=9-3', 'bytes=-0', 'bytes=0-1,3-4', 'bytes=-'):
            with self.assertRaises(ValueError): byte_range(header,100)

    def test_path_escape_and_symlink_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); allowed=root/'allowed';allowed.mkdir()
            secret=root/'other';secret.write_text('private')
            with self.assertRaises(ValueError): confined_file(secret,[allowed])
            with self.assertRaises(ValueError): confined_file(allowed/'..'/'other',[allowed])

    def test_bed_preserves_exons_and_cds_on_minus_strand(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'gene_structure.json').write_text(json.dumps({'transcripts':[
                dict(id='tx',chrom='chr1',strand='-',exons=[
                    dict(start_1based=21,end_1based=30,cds=[dict(start_1based=21,end_1based=25)]),
                    dict(start_1based=1,end_1based=10,cds=[])])]}))
            _,bed=annotations(root)
            self.assertEqual(bed.strip().split('\t'),['chr1','0','30','tx','0','-','20','25','0','2','10,10','0,20'])
