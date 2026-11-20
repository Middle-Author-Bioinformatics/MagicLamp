"""Shared behavior and real-tool checks for every maintained command; no network."""
import copy
import csv
import io as textio
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from latest import context, engine, io
from latest.models import load_models, combined_profiles
from latest.registry import GENIES, resolve_name

HMMER = all(shutil.which(x) for x in ['hmmsearch', 'hmmemit', 'hmmstat'])
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    BIO = True
except ImportError:
    BIO = False


def read_table(path, delimiter=','):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def emit(model):
    result = subprocess.run(['hmmemit', '-c', str(model.path)], capture_output=True, text=True, check=True)
    return ''.join(line for line in result.stdout.splitlines() if not line.startswith('>'))


def gbk(path, sequences, second_record=True, pseudo_gap=False):
    records = []
    for index in range(2 if second_record else 1):
        record = SeqRecord(Seq('N' * (max(1, len(sequences)) * 30000)), id='contig' + str(index), description='fixture')
        record.annotations['molecule_type'] = 'DNA'
        for rank, (identifier, sequence) in enumerate(sequences.items()):
            start = rank * 30000
            record.features.append(SeqFeature(FeatureLocation(start, start + len(sequence) * 3, strand=-1 if index else 1),
                                             type='CDS', qualifiers={'locus_tag': [identifier], 'translation': [sequence]}))
            if pseudo_gap and rank == 0:
                record.features.append(SeqFeature(FeatureLocation(start + len(sequence) * 3 + 1, start + len(sequence) * 3 + 30),
                                                 type='CDS', qualifiers={'pseudo': ['']}))
        records.append(record)
    SeqIO.write(records, str(path), 'genbank')


class Temporary(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='magiclamp test ')
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def cli(self, name, inputs, extension, mode, output='output', *extra):
        env = dict(os.environ)
        for _, (_, variable, _) in GENIES.items():
            if variable:
                env.pop(variable, None)
        result = subprocess.run([sys.executable, str(ROOT / 'MagicLamp.py'), name,
                                 '-bin_dir', str(inputs), '-bin_ext', extension, mode,
                                 '-out', str(self.root / output), *extra],
                                cwd=self.root, env=env, text=True, capture_output=True)
        return result, self.root / output


class CoreTests(Temporary):
    def test_catalogue_and_exclusions(self):
        self.assertEqual(len(GENIES), 15)
        self.assertEqual(resolve_name('iron'), 'FeGenie')
        for value in ['AbxGenie', 'abx', 'antibiotics', 'BTEXGenie', 'btex']:
            with self.assertRaisesRegex(ValueError, 'excluded'):
                resolve_name(value)
        result = subprocess.run([sys.executable, str(ROOT / 'MagicLamp.py'), 'OmniGenie', '-genie', 'btex',
                                 '-bin_dir', str(self.root), '-bin_ext', 'faa', '--orfs'], capture_output=True)
        self.assertNotEqual(result.returncode, 0)

    def test_manifest_integrity_every_library(self):
        total = 0
        for name, (directory, _, _) in GENIES.items():
            if directory:
                models, _ = load_models(name)
                self.assertEqual(len(models), len({m.model_id for m in models}))
                self.assertTrue(all(m.annotation_source and m.category for m in models))
                self.assertTrue(all(m.policy in {'ga', 'tc', 'score', 'classifier'} for m in models))
                total += len(models)
        self.assertEqual(total, 944)

    def test_override_requires_reviewed_checksums(self):
        library = self.root / 'gram'
        shutil.copytree(ROOT / 'hmms/gram', library)
        with (library / 'PF00395.hmm').open('a') as handle:
            handle.write('\n')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            load_models('GramGenie', library)

    def test_custom_threshold_is_explicit(self):
        source = (ROOT / 'hmms/gram/PF01476.hmm').read_text()
        library = self.root / 'custom';library.mkdir()
        path = library / 'nested.HMM'
        path.write_text('\n'.join(x for x in source.splitlines() if not x.startswith(('GA ', 'TC '))) + '\n')
        with self.assertRaisesRegex(ValueError, 'explicitly supply'):
            load_models('HmmGenie', library)
        models, _ = load_models('HmmGenie', library, evalue=1e-6)
        self.assertEqual(models[0].policy, 'evalue')
        with self.assertRaisesRegex(ValueError, 'No GA cutoffs'):
            load_models('HmmGenie', library, cutoff='ga')

    def test_fasta_ids_validation_and_no_assumed_positions(self):
        path = self.root / 'proteins.faa';path.write_text('>sp|123|name_99 description\nMAAA*\n')
        proteins = io.protein_fasta(path)
        self.assertEqual(proteins[0].identifier, 'sp|123|name_99')
        self.assertEqual(proteins[0].sequence, 'MAAA')
        self.assertEqual(proteins[0].order, '')
        for text, message in [('>x\nMA\n>x\nMA\n', 'Duplicate'), ('', 'No protein'), ('>x\nMA*AA\n', 'invalid')]:
            path.write_text(text)
            with self.assertRaisesRegex(ValueError, message):
                io.protein_fasta(path)

    @unittest.skipUnless(BIO, 'Biopython required')
    def test_genbank_coordinates_gaps_and_record_isolation(self):
        path = self.root / 'fixture.gbk';gbk(path, {'first': 'MAAA', 'second': 'MCCC'}, pseudo_gap=True)
        proteins = io.genbank(path)
        self.assertEqual([p.order for p in proteins], [1, 3, 1, 3])
        self.assertEqual(proteins[0].start, 1)
        self.assertEqual(proteins[0].end, 12)
        lookup = {p.identifier: p for p in proteins}
        hits = [{'target': key} for key in lookup]
        self.assertEqual(len(context.assign_clusters(hits, lookup, True, 0)), 4)
        self.assertEqual(len(context.assign_clusters(hits, lookup, True, 1)), 2)
        self.assertEqual(context.assign_clusters(hits, lookup, False), [])
        self.assertTrue(all(h['cluster_id'] == '' for h in hits))

    def test_domain_thresholds_repeat_and_multiple_models(self):
        models, _ = load_models('GramGenie')
        model = copy.copy(next(m for m in models if m.model_id == 'PF01476'))
        model.sequence_cutoff, model.domain_cutoff = 20.9, 11.8
        other = copy.copy(model);other.key = 'm999999';other.model_id = 'other'
        def line(full, domain, start, key=model.key):
            return f'p00000001 - 200 {key} - 44 1e-10 {full} 0 1 3 1e-10 1e-10 {domain} 0 1 44 {start} {start+43} {start} {start+43} 0.99 desc\n'
        def fake_search(command, log, manifest):
            Path(command[command.index('--domtblout') + 1]).write_text(
                line(20.8, 12, 1) + line(30, 11.7, 1) + line(20.9, 11.8, 1) + line(30, 12, 80) + line(30, 12, 100, other.key))
        with mock.patch.object(engine, 'checked', side_effect=fake_search):
            domains = engine.search([model, other], {('ga', 'HMMER3f'): self.root / 'models'}, self.root / 'fa',
                                    self.root, 'sample', 1, {'p00000001': io.Protein('id', 'M'*200)}, {'commands': []})
        self.assertEqual(len(domains), 3)
        hits = engine.collapse(domains)
        self.assertEqual({h['model'].model_id: h['domain_count'] for h in hits}, {'PF01476': 2, 'other': 1})

    def test_context_requires_coordinates_and_distinct_components(self):
        models, _ = load_models('FeGenie')
        selected = [m for m in models if m.path.stem in {'EetA', 'EetB', 'Ndh2', 'FmnB', 'FmnA'}]
        self.assertEqual(len(selected), 5)
        proteins = {str(i): io.Protein(str(i), 'MAAA', 'c', i*100+1, i*100+12, 1, i+1) for i in range(5)}
        hits = [{'target': str(i), 'model': m} for i, m in enumerate(selected)]
        context.assign_clusters(hits, proteins, True, 0);context.annotate('FeGenie', hits, proteins, True)
        self.assertTrue(all(h['context_status'] == 'supported_neighborhood' for h in hits))
        hits = hits[:4] + [dict(hits[0])]
        context.annotate('FeGenie', hits, proteins, True)
        self.assertTrue(all(h['context_status'] == 'insufficient_context' for h in hits))
        context.annotate('FeGenie', hits, proteins, False)
        self.assertTrue(all(h['context_status'] == 'not_assessed' for h in hits))

    def test_rnf_completeness_is_distinct_and_local(self):
        models, _ = load_models('RnfGenie')
        proteins = {str(i): io.Protein(str(i), 'MAAA', 'contig', 1+i*100, 12+i*100, 1, i+1) for i in range(6)}
        hits = [{'target': str(i), 'model': m} for i, m in enumerate(models)]
        context.assign_clusters(hits, proteins, True, 0)
        context.annotate('RnfGenie', hits, proteins, True)
        self.assertTrue(all(h['complex_completeness'] == '6/6' for h in hits))
        context.annotate('RnfGenie', hits[:3], proteins, True)
        self.assertTrue(all(h['complex_completeness'] == '3/6' for h in hits[:3]))
        context.annotate('RnfGenie', hits, proteins, False)
        self.assertTrue(all(h['complex_completeness'] == '' for h in hits))

    def test_ion_classification_missing_alternative_and_type_gate(self):
        models, _ = load_models('ATPGenie')
        model = next(m for m in models if m.accession.split('.')[0] == 'TIGR01260')
        fmodel = next(m for m in models if m.name == 'F-type_H_pumping')
        vmodel = next(m for m in models if m.name == 'V-type_Na_binding')
        hits = [{'target': 'candidate', 'model': model}]
        lookup = {'candidate': io.Protein('atpE', 'MAAA'), 'unrelated': io.Protein('unrelated', 'MAAA')}
        args = engine.make_parser('ATPGenie').parse_args(['-bin_dir', str(self.root), '-bin_ext', 'faa', '--orfs'])
        with mock.patch.object(engine, 'search', return_value=[{'target': 'candidate', 'model': fmodel, 'bit_score': 20},
                                                              {'target': 'candidate', 'model': vmodel, 'bit_score': 100}]):
            engine.ion_preferences(hits, models, {}, self.root/'all.faa', self.root, 'test', args, lookup, {'commands': []})
        self.assertEqual(hits[0]['ion_model'], 'F-type_H_pumping')
        self.assertEqual(hits[0]['ion_confidence'], 'unresolved')
        self.assertEqual(hits[0]['ion_score_delta'], '')
        self.assertNotIn('unrelated', (self.root/'test.proteolipids.faa').read_text())

    def test_failed_search_manifest_and_no_false_success(self):
        inputs = self.root / 'inputs';inputs.mkdir();(inputs/'a.faa').write_text('>a\nMAAA\n')
        out = self.root / 'failed'
        with mock.patch.object(engine.shutil, 'which', return_value='hmmsearch'), \
             mock.patch.object(engine, 'tool_version', return_value='test'), \
             mock.patch.object(engine, 'checked', side_effect=subprocess.CalledProcessError(1, ['hmmsearch'])):
            self.assertEqual(engine.main('GramGenie', ['-bin_dir', str(inputs), '-bin_ext', 'faa', '--orfs', '-out', str(out)]), 1)
        self.assertEqual(json.loads((out/'run.json').read_text())['status'], 'failed')
        self.assertFalse((out/'summary.csv').exists())
        self.assertEqual(engine.main('GramGenie', ['-bin_dir', str(inputs), '-bin_ext', 'faa', '--orfs', '-out', str(out)]), 1)

    def test_strict_flags(self):
        inputs = self.root / 'input';inputs.mkdir();(inputs/'a.faa').write_text('>a\nMAAA\n')
        for extra in [('--orfs', '--gbk'), ('--not-real',), ('--meta', '--orfs')]:
            result = subprocess.run([sys.executable, str(ROOT/'MagicLamp.py'), 'GramGenie', '-bin_dir', str(inputs),
                                     '-bin_ext', 'faa', '-out', str(self.root/'out'), *extra], capture_output=True)
            self.assertNotEqual(result.returncode, 0)

    def test_setup_is_cwd_independent_and_non_destructive(self):
        result = subprocess.run(['bash', str(ROOT/'setup.sh'), '--dry-run', '--name', 'test-env'],
                                cwd=self.root, capture_output=True, text=True, check=True)
        self.assertIn(str(ROOT/'environment.yml'), result.stdout)
        prefix = self.root/'fake env';(prefix/'conda-meta').mkdir(parents=True)
        bindir = self.root/'bin';bindir.mkdir()
        conda = bindir/'conda'
        conda.write_text("""#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$FAKE_CALLS"
if [[ "$1" == env ]]; then exit 17; fi
if [[ "$*" == *'print(sys.prefix)'* ]]; then printf '%s\\n' "$FAKE_PREFIX"; fi
""")
        conda.chmod(0o755)
        calls = self.root/'calls.txt'
        env = dict(os.environ, PATH=str(bindir)+':'+os.environ['PATH'], FAKE_PREFIX=str(prefix), FAKE_CALLS=str(calls))
        failed = subprocess.run(['bash', str(ROOT/'setup.sh')], env=env, capture_output=True)
        self.assertEqual(failed.returncode, 17)
        creation = calls.read_text()
        self.assertIn('env create --name magiclamp --file ', creation)
        self.assertNotIn('--yes', creation)
        calls.write_text('')
        subprocess.run(['bash', str(ROOT/'setup.sh'), '--install-only'], env=env, check=True, capture_output=True)
        recorded = calls.read_text()
        self.assertIn('python -m pip install --no-deps --no-build-isolation', recorded)
        self.assertIn('magiclamp doctor --smoke-test', recorded)
        self.assertNotIn('env create', recorded)
        self.assertFalse((prefix/'etc/conda/activate.d').exists())
        (prefix/'bin').mkdir()
        launcher = prefix/'bin/MagicLamp.py'
        launcher.write_text('#!/bin/sh\n# unrelated executable\n')
        refused = subprocess.run(['bash', str(ROOT/'setup.sh'), '--install-only'], env=env, capture_output=True)
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn('unrelated executable', launcher.read_text())

    def test_setup_supports_legacy_conda_env_parser(self):
        prefix = self.root/'environment'
        (prefix/'conda-meta').mkdir(parents=True)
        bindir = self.root/'bin';bindir.mkdir()
        conda = bindir/'conda'
        conda.write_text("""#!/usr/bin/env bash
if [[ "$1" == env && "$2" == create ]]; then
    for argument in "$@"; do
        if [[ "$argument" == --yes || "$argument" == -y ]]; then
            echo 'conda-env: error: unrecognized arguments: --yes' >&2
            exit 2
        fi
    done
    [[ "$CONDA_ALWAYS_YES" == true ]] || exit 18
elif [[ "$*" == *'print(sys.prefix)'* ]]; then
    printf '%s\\n' "$FAKE_PREFIX"
else
    # The setting must apply only to env creation, not leak into later calls.
    [[ ! ${CONDA_ALWAYS_YES+x} ]] || exit 19
fi
""")
        conda.chmod(0o755)
        env = dict(os.environ, PATH=str(bindir)+':'+os.environ['PATH'], FAKE_PREFIX=str(prefix))
        env.pop('CONDA_ALWAYS_YES', None)
        result = subprocess.run(['bash', str(ROOT/'setup.sh'), '--plots'], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)



@unittest.skipUnless(HMMER, 'HMMER executables required')
class IntegrationTests(Temporary):
    def test_every_runtime_profile_is_parseable(self):
        for name, (directory, _, _) in GENIES.items():
            if directory:
                models, _ = load_models(name)
                dest = self.root/name;dest.mkdir()
                for path in combined_profiles(models, dest).values():
                    subprocess.run(['hmmstat', str(path)], check=True, capture_output=True)

    @unittest.skipUnless(BIO and shutil.which('blastp'), 'Biopython and BLASTP required')
    def test_all_commands_proteins_and_genbank_with_zero_sample(self):
        for name, (directory, _, _) in GENIES.items():
            effective = 'GramGenie' if name in {'HmmGenie', 'OmniGenie'} else name
            models, _ = load_models(effective)
            chosen = next((m for m in models if m.accession.split('.')[0] == 'TIGR01260'), models[0])
            sequence = emit(chosen)
            extras = ['-genie', 'gram'] if name == 'OmniGenie' else ['-hmm_dir', str(ROOT/'hmms/gram')] if name == 'HmmGenie' else []
            for mode in ['--orfs', '--gbk']:
                with self.subTest(genie=name, mode=mode):
                    inputs = self.root/(name+mode);inputs.mkdir()
                    ext = 'faa' if mode == '--orfs' else 'gbk'
                    if ext == 'faa':
                        (inputs/f'positive, sample.{ext}').write_text('>sp|positive_99|id description\n'+sequence+'\n')
                        (inputs/f'zero.{ext}').write_text('>negative_1\n'+'A'*150+'\n')
                    else:
                        gbk(inputs/f'positive, sample.{ext}', {'positive': sequence})
                        gbk(inputs/f'zero.{ext}', {'negative': 'A'*150}, second_record=False)
                    result, out = self.cli(name, inputs, ext, mode, name+mode+'out', *extras)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    manifest = json.loads((out/'run.json').read_text())
                    self.assertEqual(manifest['status'], 'complete')
                    self.assertEqual(manifest['clustering'], mode == '--gbk')
                    self.assertFalse(any(c[0] == 'prodigal' for c in manifest['commands']))
                    rows = read_table(out/'summary.csv')
                    self.assertTrue(rows, f'No positive control hit: {name}')
                    self.assertTrue(any(r['model'] == chosen.model_id for r in rows))
                    self.assertTrue(all(bool(r['cluster_id']) == (mode == '--gbk') for r in rows))
                    if mode == '--orfs':
                        self.assertTrue(all(r['protein_id'] == 'sp|positive_99|id' and not r['contig'] for r in rows))
                    else:
                        self.assertEqual({r['contig'] for r in rows}, {'contig0','contig1'})
                        self.assertEqual(len(read_table(out/'clusters.tsv', '\t')), 2)
                    with (out/'counts.csv').open() as handle:
                        matrix = list(csv.reader(handle))
                    column = matrix[0].index('zero.'+ext)
                    self.assertTrue(all(row[column] == '0' for row in matrix[1:]))
                    self.assertEqual(len(read_table(out/'samples.tsv', '\t')), 2)
                    if name == 'FeGenie':
                        self.assertEqual(len([c for c in manifest['commands'] if c[0] == 'blastp']), 4)
                        self.assertEqual(len(manifest['references']), 2)
                    if name == 'ATPGenie':
                        candidates = [r for r in rows if r['accession'].split('.')[0] == 'TIGR01260']
                        self.assertTrue(all(r['ion_model'].startswith('F-type') for r in candidates))

    def test_gram_repeats_multiple_families_and_normalization(self):
        models, _ = load_models('GramGenie')
        inputs = self.root/'proteins';inputs.mkdir()
        sequences = {m.model_id: (emit(m)+'GGGGSGGGGS')*3 for m in models}
        sequences['multi_domain'] = sequences['PF01476'] + sequences['PF07501']
        (inputs/'positive.faa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in sequences.items()))
        (inputs/'zero.faa').write_text('>negative\n'+'A'*150+'\n')
        result, out = self.cli('GramGenie', inputs, 'faa', '--orfs', 'repeat out', '--norm')
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = read_table(out/'summary.csv')
        for model in models:
            hit = next(r for r in rows if r['protein_id'] == model.model_id and r['model'] == model.model_id)
            self.assertGreaterEqual(int(hit['domain_count']), 2)
        self.assertTrue({'PF01476','PF07501'} <= {r['model'] for r in rows if r['protein_id'] == 'multi_domain'})
        with (out/'counts.csv').open() as handle: matrix=list(csv.reader(handle))
        for row in matrix[1:]:
            count=sum(r['model'] == row[0] for r in rows)
            self.assertAlmostEqual(float(row[1]), count*1000/len(sequences))

    @unittest.skipUnless(BIO, 'Biopython required')
    def test_genbank_no_cluster_override(self):
        inputs=self.root/'gbk';inputs.mkdir()
        models,_=load_models('GramGenie');gbk(inputs/'sample.gbk', {'positive': emit(models[0])})
        result,out=self.cli('GramGenie', inputs, 'gbk', '--gbk', 'nocluster', '--no-cluster')
        self.assertEqual(result.returncode, 0, result.stderr)
        rows=read_table(out/'summary.csv')
        self.assertTrue(rows)
        self.assertTrue(all(r['contig'] and not r['cluster_id'] and r['context_status']=='not_assessed' for r in rows))
        self.assertEqual(read_table(out/'clusters.tsv','\t'), [])

    @unittest.skipUnless(__import__('importlib.util').util.find_spec('matplotlib'), 'matplotlib required')
    def test_plot_uses_count_matrix(self):
        inputs=self.root/'plot';inputs.mkdir();(inputs/'zero.faa').write_text('>a\n'+'A'*150+'\n')
        result,out=self.cli('GramGenie', inputs, 'faa', '--orfs', 'plotout', '--makeplots', '--norm')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((out/'heatmap.png').read_bytes().startswith(b'\x89PNG'))

    @unittest.skipUnless(shutil.which('prodigal'), 'Prodigal required')
    def test_zero_proteins_normalization(self):
        inputs=self.root/'dna';inputs.mkdir();(inputs/'empty.fna').write_text('>contig\n'+'N'*60+'\n')
        result = subprocess.run([sys.executable, str(ROOT/'MagicLamp.py'), 'GramGenie', '-bin_dir', str(inputs),
                                 '-bin_ext', 'fna', '--meta', '--norm', '-out', str(self.root/'zeroout')], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(read_table(self.root/'zeroout/samples.tsv','\t')[0]['protein_count'], '0')
        with (self.root/'zeroout/counts.csv').open() as handle:
            self.assertTrue(all(r[1] == '' for r in list(csv.reader(handle))[1:]))


if __name__ == '__main__':
    unittest.main()
