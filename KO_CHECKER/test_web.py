"""Offline regression tests for data reuse, accession handling and web jobs."""
import csv
import io
import json
from pathlib import Path
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pipeline
from pipeline import Pipeline, atomic_json, mark_complete, valid_cache, validate_request
from webapp import Jobs, handler


class InputTests(unittest.TestCase):
    def test_optional_threads_accepts_only_bounded_integers(self):
        payload={'gene':'Phf3','wt':'SRX1','ko':'SRX2'}
        for value in (None, '', '  '):
            self.assertNotIn('threads', validate_request({**payload,'threads':value}))
        for value in (1,256,'8',' 4 '):
            self.assertEqual(validate_request({**payload,'threads':value})['threads'],int(value))
        for value in (0,-1,257,1.5,True,False,'2.5','1e2',[],{}):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError,'スレッド'):
                validate_request({**payload,'threads':value})

    def test_normalizes_multiple_experiments(self):
        request = validate_request({'gene': 'Phf3', 'wt': 'srx1, ERX2\nDRX3', 'ko': 'SRX4'})
        self.assertEqual(request['wt'], ['SRX1', 'ERX2', 'DRX3'])


    def test_rejects_overlap_duplicates_invalid_ids_and_injection(self):
        for wt, ko in [('SRX1', 'SRX1'), ('SRX1,SRX1', 'SRX2'), ('ERR123', 'SRX2'), ('SRR123', 'SRX2'), ('DRR123', 'SRX2'), ('SRS123', 'SRX2'), ('SRX1; rm -rf /', 'SRX2'), ('', 'SRX2')]:
            with self.subTest(wt=wt, ko=ko), self.assertRaises(ValueError):
                validate_request({'gene': 'Phf3', 'wt': wt, 'ko': ko})


class CacheAndMetadataTests(unittest.TestCase):
    @staticmethod
    def metadata_response():
        row = {'run_accession':'SRR1', 'experiment_accession':'SRX1', 'tax_id':'10090',
               'library_strategy':'RNA-Seq', 'library_layout':'SINGLE'}
        return io.BytesIO(json.dumps([row]).encode())

    def test_transient_ena_error_is_retried(self):
        error = urllib.error.HTTPError('https://ena.example', 500, 'Internal Server Error', {}, None)
        with tempfile.TemporaryDirectory() as temp, \
             patch('pipeline.urllib.request.urlopen',side_effect=[error,self.metadata_response()]) as fetch, \
             patch('pipeline.time.sleep') as sleep:
            rows=pipeline.resolve_runs('SRX1','mm10',temp)
            self.assertEqual(rows[0]['run_accession'],'SRR1')
            self.assertEqual(fetch.call_count,2)
            sleep.assert_called_once_with(2)

    def test_search_fallback_recovers_from_filereport_http500(self):
        error = urllib.error.HTTPError('https://ena.example',500,'Internal Server Error',{},None)
        notes=[]
        with tempfile.TemporaryDirectory() as temp, \
             patch('pipeline.urllib.request.urlopen',side_effect=[error,error,error,self.metadata_response()]) as fetch, \
             patch('pipeline.time.sleep'):
            rows=pipeline.resolve_runs('SRX1','mm10',temp,notes.append)
            self.assertEqual(rows[0]['run_accession'],'SRR1')
            self.assertEqual(fetch.call_count,4)
            self.assertIn('/search?',fetch.call_args.args[0].full_url)
            self.assertTrue(any('検索API' in n for n in notes))
            self.assertTrue((Path(temp)/'metadata'/'SRX1.json').is_file())

    def test_outage_has_actionable_error_and_does_not_cache_failure(self):
        error=urllib.error.HTTPError('https://ena.example',503,'Service Unavailable',{},None)
        with tempfile.TemporaryDirectory() as temp, \
             patch('pipeline.urllib.request.urlopen',side_effect=error) as fetch, \
             patch('pipeline.time.sleep'), self.assertRaisesRegex(RuntimeError,'SRX1.*ENA'):
            try:
                pipeline.resolve_runs('SRX1','mm10',temp)
            finally:
                self.assertEqual(fetch.call_count,6)
                self.assertFalse((Path(temp)/'metadata'/'SRX1.json').exists())

    def test_empty_report_uses_search_and_genuinely_missing_accession_is_distinct(self):
        with patch('pipeline.urllib.request.urlopen',side_effect=[io.BytesIO(b'[]'),self.metadata_response()]):
            self.assertEqual(pipeline.fetch_run_metadata('SRX1')[0]['run_accession'],'SRR1')
        with patch('pipeline.urllib.request.urlopen',side_effect=[io.BytesIO(b'[]'),io.BytesIO(b'[]')]), \
             self.assertRaisesRegex(ValueError,'公開runが見つかりません'):
            pipeline.fetch_run_metadata('SRX1')

    def test_invalid_json_response_recovers_on_retry(self):
        with patch('pipeline.urllib.request.urlopen',side_effect=[io.BytesIO(b'<html>service error</html>'),self.metadata_response()]), \
             patch('pipeline.time.sleep'):
            self.assertEqual(pipeline.fetch_run_metadata('SRX1')[0]['run_accession'],'SRR1')

    def test_partial_changed_or_wrong_reference_cache_is_not_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'run.bam'
            path.write_bytes(b'complete')
            signature = {'reference': 'mm10'}
            self.assertFalse(valid_cache(temp, signature))
            mark_complete(temp, signature, [path])
            self.assertTrue(valid_cache(temp, signature))
            self.assertFalse(valid_cache(temp, {'reference': 'mm39'}))
            path.write_bytes(b'truncated')
            self.assertFalse(valid_cache(temp, signature))

    def test_resolves_all_runs_and_uses_metadata_cache(self):
        rows = [{'run_accession': run, 'experiment_accession': 'SRX1', 'tax_id': '10090',
                 'scientific_name': 'Mus musculus', 'library_strategy': 'RNA-Seq', 'library_layout': 'PAIRED'}
                for run in ['SRR2', 'SRR1']]
        with tempfile.TemporaryDirectory() as temp, patch('pipeline.urllib.request.urlopen', return_value=io.BytesIO(json.dumps(rows).encode())) as fetch:
            self.assertEqual([r['run_accession'] for r in pipeline.resolve_runs('SRX1', 'mm10', temp)], ['SRR1', 'SRR2'])
            pipeline.resolve_runs('SRX1', 'mm10', temp)
            self.assertEqual(fetch.call_count, 1)
            with self.assertRaisesRegex(ValueError, '生物種'):
                pipeline.resolve_runs('SRX1', 'hg38', temp)

    def test_refuses_non_rna_seq(self):
        row = {'run_accession': 'SRR1', 'experiment_accession': 'SRX1', 'tax_id': '10090',
               'library_strategy': 'WGS', 'library_layout': 'PAIRED'}
        with tempfile.TemporaryDirectory() as temp, patch('pipeline.urllib.request.urlopen', return_value=io.BytesIO(json.dumps([row]).encode())):
            with self.assertRaisesRegex(ValueError, 'RNA-Seq'):
                pipeline.resolve_runs('SRX1', 'mm10', temp)




class PipelineTests(unittest.TestCase):
    def make_pipeline(self, temp):
        root = Path(temp)
        return Pipeline(root / 'data', root / 'references', 2, lambda **kw: None, root / 'pipeline.log')

    def test_fastq_download_and_reuse_and_incomplete_retry(self):
        with tempfile.TemporaryDirectory() as temp:
            engine = self.make_pipeline(temp)
            row = {'run_accession': 'SRR1', 'library_layout': 'PAIRED'}
            calls = []
            def command(args, stdout=None):
                calls.append(args)
                if args[0] == 'prefetch':
                    archive = Path(args[-1]); archive.mkdir(parents=True, exist_ok=True)
                    (archive / 'SRR1.sra').write_bytes(b'sra')
                if args[0] == 'fasterq-dump':
                    directory = Path(args[args.index('-O')+1])
                    for suffix in ('_1', '_2'):
                        (directory / ('SRR1'+suffix+'.fastq')).write_text('@r\nAC\n+\nII\n')
            with patch.object(engine, 'command', side_effect=command):
                paths = engine.fastqs(row, {'fasterq-dump': '1'})
                self.assertEqual(len(paths), 2)
                self.assertEqual(len(calls), 3)
                engine.fastqs(row, {'fasterq-dump': '1'})
                self.assertEqual(len(calls), 3)
                paths[1].write_bytes(b'')
                engine.fastqs(row, {'fasterq-dump': '1'})
                self.assertEqual(len(calls), 6)
                self.assertEqual(Path(calls[0][-1]).name, 'SRR1')

    def test_incomplete_fastq_never_gets_complete_marker(self):
        with tempfile.TemporaryDirectory() as temp:
            engine = self.make_pipeline(temp)
            def command(args, stdout=None):
                if args[0] == 'prefetch':
                    archive = Path(args[-1]); archive.mkdir(parents=True, exist_ok=True)
                    (archive / 'SRR1.sra').write_bytes(b'sra')
            with patch.object(engine, 'command', side_effect=command), self.assertRaisesRegex(RuntimeError, 'FASTQ'):
                engine.fastqs({'run_accession': 'SRR1', 'library_layout': 'PAIRED'}, {'fasterq-dump': '1'})
            self.assertFalse((engine.data / 'fastq' / 'SRR1' / 'complete.json').exists())

    def test_fastq_accepts_direct_and_nested_prefetch_archives(self):
        for nested in (False, True):
            for extension in ('.sra', '.sralite'):
                with self.subTest(nested=nested, extension=extension), tempfile.TemporaryDirectory() as temp:
                    engine = self.make_pipeline(temp)
                    destination = engine.data / 'sra' / 'SRR1'
                    if nested:
                        destination = destination / 'SRR1'
                    archive = destination / ('SRR1' + extension)
                    def command(args, stdout=None):
                        if args[0] == 'prefetch':
                            archive.parent.mkdir(parents=True, exist_ok=True)
                            archive.write_bytes(b'sra')
                        elif args[0] == 'vdb-validate':
                            self.assertEqual(args[1], archive)
                        elif args[0] == 'fasterq-dump':
                            self.assertEqual(args[1], archive)
                            directory = Path(args[args.index('-O') + 1])
                            (directory / 'SRR1.fastq').write_text('@r\nAC\n+\nII\n')
                    with patch.object(engine, 'command', side_effect=command):
                        result = engine.fastqs({'run_accession':'SRR1','library_layout':'SINGLE'}, {'fasterq-dump':'v1'})
                    self.assertEqual(len(result), 1)
                    self.assertTrue((result[0].parent / 'complete.json').exists())

    def test_completed_bam_skips_fastq_and_mapping(self):
        with tempfile.TemporaryDirectory() as temp:
            engine = self.make_pipeline(temp)
            index = Path(temp) / 'index' / 'genome'
            index.parent.mkdir()
            atomic_json(index.parent / 'complete.json', {'reference': 'mm10'})
            row = {'run_accession': 'SRR1', 'library_layout': 'SINGLE'}
            versions = {'hisat2': 'v1', 'samtools': 'v2'}
            signature = {'version': pipeline.PIPELINE_VERSION, 'run': 'SRR1', 'layout': 'SINGLE',
                         'index': {'reference': 'mm10'}, 'hisat2': 'v1', 'samtools': 'v2',
                         'options': ['--dta', '--rg-id', 'SRR1']}
            directory = engine.data / 'bam' / 'SRR1' / pipeline.digest(signature)
            directory.mkdir(parents=True)
            bam = directory / 'SRR1.bam'; bam.write_bytes(b'bam')
            bai = Path(str(bam) + '.bai'); bai.write_bytes(b'index')
            mark_complete(directory, signature, [bam, bai])
            with patch.object(engine, 'command') as command, patch.object(engine, 'fastqs') as fastqs:
                self.assertEqual(engine.mapped_run(row, index, versions), bam)
                fastqs.assert_not_called()
                command.assert_called_once_with(['samtools', 'quickcheck', bam])

    def test_alignment_failure_does_not_publish_bam(self):
        with tempfile.TemporaryDirectory() as temp:
            engine = self.make_pipeline(temp)
            index = Path(temp) / 'index' / 'genome'; index.parent.mkdir()
            atomic_json(index.parent / 'complete.json', {'reference': 'mm10'})
            from unittest.mock import MagicMock
            aligner = MagicMock(); aligner.wait.return_value = 1; aligner.poll.return_value = 1
            sorter = MagicMock(); sorter.wait.return_value = 0
            with patch.object(engine, 'fastqs', return_value=[Path(temp)/'read.fastq']), \
                 patch('pipeline.subprocess.Popen', side_effect=[aligner,sorter]), \
                 self.assertRaisesRegex(RuntimeError, 'HISAT2=1'):
                engine.mapped_run({'run_accession':'SRR1','library_layout':'SINGLE'},index,{'hisat2':'v1','samtools':'v2'})
            self.assertEqual(list(engine.data.rglob('complete.json')), [])

    def test_observation_export_distinguishes_zero_alt_and_no_depth(self):
        from unittest.mock import MagicMock
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'events.csv').write_text('event_id,chrom,pos_1based,ref,alt\ne1,chr1,10,A,AT\n')
            atomic_json(root/'run.json',{'thresholds':{'min_mapq':0,'min_bq':20,'min_event_depth':10}})
            samples=[{'srx_id':sid,'metadata_label':'WT','bam_path':'unused'} for sid in ['SRX1','SRX2']]
            fake=SimpleNamespace(FastaFile=MagicMock(),AlignmentFile=MagicMock())
            tx=SimpleNamespace(chrom='chr1')
            with patch('pipeline.rnako.pysam',fake), patch('pipeline.rnako.load_gene',return_value={'t1':tx}), \
                 patch('pipeline.rnako.resolve_contig',return_value='chr1'), \
                 patch('pipeline.rnako.observe_event',side_effect=[(100,0),(0,0)]):
                pipeline.write_observations(samples,root,'unused','unused','Phf3')
            with (root/'observations.csv').open() as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(rows[0]['local_alt_fraction'],'0.0')
            self.assertEqual(rows[0]['coverage_status'],'covered')
            self.assertEqual(rows[1]['local_alt_fraction'],'')
            self.assertEqual(rows[1]['coverage_status'],'low_depth')

    def test_pipeline_orchestrates_multiple_runs_without_changing_core_arguments(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); engine = self.make_pipeline(temp)
            job = root / 'job'; job.mkdir()
            fasta = root / 'reference.fa'; fasta.write_text('>chr1\nATG\n')
            gtf = root / 'reference.gtf'; gtf.write_text('annotation')
            def runs(accession, assembly, cache, note=None):
                return [{'run_accession': 'SRR'+accession[3:]+str(i), 'experiment_accession':accession} for i in range(2 if accession=='SRX1' else 1)]
            def mapped(row, index, versions):
                bam=root/(row['run_accession']+'.bam');bam.write_bytes(b'bam');return bam
            states=[];engine.update=lambda **kw: states.append(kw)
            with patch.object(engine, 'tools', return_value={'hisat2':'version'}), \
                 patch('pipeline.resolve_runs', side_effect=runs), \
                 patch('pipeline.rnako.prepare_reference', return_value=(fasta,gtf)), \
                 patch('pipeline.rnako.load_gene'), patch.object(engine,'index',return_value=root/'index'), \
                 patch.object(engine,'mapped_run',side_effect=mapped), \
                 patch.object(engine,'sample_bam',side_effect=lambda a,b,v:b[0]) as aggregate, \
                 patch.object(engine,'command') as command, patch('pipeline.write_observations'), \
                 patch('pipeline.write_exon_coverage') as coverage:
                engine.run(validate_request({'gene':'Phf3','wt':'SRX1','ko':'SRX2'}),job)
            self.assertEqual(len(aggregate.call_args_list[0].args[1]), 2)
            self.assertEqual(len(aggregate.call_args_list[1].args[1]), 1)
            args=command.call_args.args[0]
            self.assertIn('--samples',args);self.assertIn('Phf3',args)
            with (job/'input_samples.csv').open() as handle:
                samples=list(csv.DictReader(handle))
            self.assertEqual([s['metadata_label'] for s in samples],['WT','KO'])
            self.assertEqual(states[-1]['status'],'completed')
            coverage.assert_called_once()
            self.assertEqual(json.loads((job/'provenance.json').read_text())['threads'],engine.threads)
            self.assertEqual(len(json.loads((job/'provenance.json').read_text())['experiments']),2)




class WebTests(unittest.TestCase):
    def test_http_submit_poll_results_download_and_origin_checks(self):
        with tempfile.TemporaryDirectory() as temp:
            def run(engine, request, directory):
                self.assertEqual(engine.threads, 2)
                output=directory/'output';output.mkdir()
                for name in ['samples','events','junctions','consequences','observations','exon_coverage']:
                    (output/(name+'.csv')).write_text('srx_id,predicted_class\nSRX1,WT_like\n',encoding='utf-8')
                engine.update(status='completed',stage='completed',progress=100,message='done')
            with patch.object(Pipeline,'run',run):
                jobs=Jobs(Path(temp)/'jobs',Path(temp)/'data',Path(temp)/'refs',1)
                server=ThreadingHTTPServer(('127.0.0.1',0),handler(jobs))
                thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
                base=f'http://127.0.0.1:{server.server_port}'
                try:
                    with urllib.request.urlopen(base+'/') as response:
                        self.assertIn('KO_Checker',response.read().decode())
                    payload=json.dumps({'gene':'Phf3','wt':'SRX1','ko':'SRX2','threads':'2'}).encode()
                    request=urllib.request.Request(base+'/api/jobs',data=payload,headers={'Content-Type':'application/json'})
                    with urllib.request.urlopen(request) as response:
                        self.assertEqual(response.status,202); job=json.load(response)
                        self.assertEqual(job['threads'],2)
                    jobs.queue.join()
                    with urllib.request.urlopen(base+'/api/jobs/'+job['id']) as response:
                        self.assertEqual(json.load(response)['status'],'completed')
                    with urllib.request.urlopen(base+'/api/jobs/'+job['id']+'/results') as response:
                        result=json.load(response)
                        self.assertEqual(result['tables']['samples'][0]['predicted_class'],'WT_like')
                        self.assertEqual(result['counts']['exon_coverage'],1)
                    with urllib.request.urlopen(base+'/api/jobs/'+job['id']+'/files/exon_coverage.csv') as response:
                        self.assertIn(b'SRX1',response.read())
                    (Path(temp)/'jobs'/job['id']/'output'/'exon_coverage.csv').unlink()
                    with urllib.request.urlopen(base+'/api/jobs/'+job['id']+'/results') as response:
                        self.assertEqual(json.load(response)['tables']['exon_coverage'],[])
                    with urllib.request.urlopen(base+'/api/jobs/'+job['id']+'/files/samples.csv') as response:
                        self.assertIn(b'SRX1',response.read())
                    request.add_header('Origin','https://unrelated.example')
                    with self.assertRaises(urllib.error.HTTPError) as caught:
                        urllib.request.urlopen(request)
                    self.assertEqual(caught.exception.code,403)
                    with self.assertRaises(urllib.error.HTTPError) as caught:
                        urllib.request.urlopen(base+'/api/jobs/'+job['id']+'/files/..%2Fjob.json')
                    self.assertEqual(caught.exception.code,404)
                finally:
                    server.shutdown();server.server_close();thread.join()

    def test_restart_marks_running_jobs_interrupted(self):
        with tempfile.TemporaryDirectory() as temp:
            jid='a'*32
            atomic_json(Path(temp)/jid/'job.json',{'id':jid,'status':'running','created_at':'2026-01-01'})
            jobs=Jobs(temp,'unused','unused',1)
            self.assertEqual(jobs.get(jid)['status'],'paused')
            self.assertEqual(jobs.get(jid)['stage'],'interrupted')


if __name__ == '__main__':
    unittest.main()
