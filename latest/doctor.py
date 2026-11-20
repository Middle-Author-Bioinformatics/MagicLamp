"""Offline dependency/model checks and an end-to-end GramGenie installation test."""
import argparse
import csv
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from . import __version__
from .models import load_models
from .registry import GENIES, HMM_ROOT

TOOLS = {
    'hmmsearch': (['-h'], (3, 3)),
    'prodigal': (['-v'], (2, 6, 3)),
    'blastp': (['-version'], (2, 12)),
}


def check_tool(executable, flags, minimum):
    path = shutil.which(executable)
    if not path:
        raise ValueError('{} missing from PATH'.format(executable))
    result = subprocess.run([path, *flags], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    output = result.stdout
    patterns = {'hmmsearch': r'HMMER\s+(\d+(?:\.\d+)+)', 'prodigal': r'Prodigal\s+V?(\d+(?:\.\d+)+)',
                'blastp': r'blastp:\s*(\d+(?:\.\d+)+)'}
    match = re.search(patterns[executable], output, re.I)
    if not match:
        raise ValueError('Cannot identify {} version from {}'.format(executable, path))
    version = tuple(int(x) for x in match[1].split('.'))
    padded = version + (0,) * max(0, len(minimum) - len(version))
    if padded < minimum:
        raise ValueError('{} {} is too old; need {}'.format(executable, match[1], '.'.join(map(str, minimum))))
    return '{} ({})'.format(match[1], path)


def smoke_test():
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from .engine import main as run
    models, _ = load_models('GramGenie')
    if not shutil.which('hmmemit'):
        raise ValueError('hmmemit is needed for the installation smoke test; install the full HMMER package')
    with tempfile.TemporaryDirectory(prefix='magiclamp installation ') as directory:
        root = Path(directory)
        proteins = root/'proteins';proteins.mkdir()
        sequences = {}
        for model in models:
            emitted = subprocess.run(['hmmemit', '-c', str(model.path)], text=True, capture_output=True, check=True)
            sequences[model.model_id] = ''.join(line for line in emitted.stdout.splitlines() if not line.startswith('>'))
        (proteins/'positive.faa').write_text(''.join('>{}\n{}\n'.format(k, v) for k, v in sequences.items()))
        (proteins/'zero.faa').write_text('>zero\n' + 'A'*150 + '\n')
        dna = root/'genomes';dna.mkdir()
        shutil.copyfile(Path(__file__).parent/'fixtures/B_subtilis_first_30kb.fna', dna/'fixture.fna')
        gbk = root/'genbank';gbk.mkdir()
        sequence = sequences['PF01476']
        record = SeqRecord(Seq('N'*(len(sequence)*3)), id='test_contig', description='installation fixture')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(FeatureLocation(0, len(record), strand=1), type='CDS',
                                      qualifiers={'locus_tag':['lysM_control'], 'translation':[sequence]})]
        SeqIO.write([record], str(gbk/'fixture.gbk'), 'genbank')
        for source, extension, flags in [(dna, 'fna', ['--meta']), (proteins, 'faa', ['--orfs']), (gbk, 'gbk', ['--gbk'])]:
            out = root/(extension+'_results')
            if run('GramGenie', ['-bin_dir', str(source), '-bin_ext', extension, '-out', str(out), *flags]):
                raise ValueError('GramGenie {} smoke test failed'.format(extension))
            manifest = json.loads((out/'run.json').read_text())
            with (out/'summary.csv').open() as handle:
                hits = list(csv.DictReader(handle))
            if manifest['status'] != 'complete' or not hits:
                raise ValueError('Incomplete/no-hit {} positive control'.format(extension))
            if extension == 'fna' and manifest['samples'][0]['protein_count'] != 25:
                raise ValueError('Nucleotide fixture should produce 25 proteins')
            if extension == 'faa':
                own = {r['model'] for r in hits if r['protein_id'] == r['model']}
                if own != {m.model_id for m in models} or any(r['cluster_id'] for r in hits):
                    raise ValueError('Protein controls did not retain all six models without clustering')
                with (out/'counts.csv').open() as handle:
                    matrix = list(csv.reader(handle))
                zero = matrix[0].index('zero.faa')
                if not all(row[zero] == '0' for row in matrix[1:]):
                    raise ValueError('Zero-hit sample was not preserved')
            elif not all(row['cluster_id'] for row in hits):
                raise ValueError('Coordinate-supported input has unclustered hits')
            print('PASS GramGenie {} input'.format(extension))


def main(argv=None):
    parser = argparse.ArgumentParser(prog='magiclamp doctor', description=__doc__)
    parser.add_argument('--plots', action='store_true', help='Also check optional plotting dependencies')
    parser.add_argument('--smoke-test', action='store_true', help='Run offline positive controls for all three GramGenie input modes')
    args = parser.parse_args(argv)
    failures = []
    def check(label, action):
        try:
            detail = action()
            print('PASS {}{}'.format(label, ': '+str(detail) if detail is not None else ''))
        except (OSError, ValueError, ImportError, importlib.metadata.PackageNotFoundError, subprocess.CalledProcessError) as exc:
            failures.append(label)
            print('FAIL {}: {}'.format(label, exc), file=sys.stderr)
    print('MagicLamp {} installation check'.format(__version__))
    def python_version():
        if not (3,11) <= sys.version_info[:2] < (3,14):
            raise ValueError('Python 3.11–3.13 required')
        return sys.version.split()[0]
    check('Python', python_version)
    for executable, (flags, minimum) in TOOLS.items():
        check(executable, lambda e=executable, f=flags, m=minimum: check_tool(e, f, m))
    for module, distribution in [('Bio','biopython')] + ([('numpy','numpy'), ('matplotlib','matplotlib')] if args.plots else []):
        def package(m=module, d=distribution):
            importlib.import_module(m)
            return importlib.metadata.version(d)
        check(distribution, package)
    print('Bundled model location: {}'.format(HMM_ROOT))
    for name, (directory, _, _) in GENIES.items():
        if directory:
            check(name, lambda n=name: '{} validated profiles'.format(len(load_models(n)[0])))
    for name in ['TherJR_SLCs.faa', 'geobacter_PCCs.faa']:
        def reference(n=name):
            library = Path(os.environ.get(GENIES['FeGenie'][1]) or HMM_ROOT/'iron')
            path = library/'iron_reduction/non-aligned'/n
            if not path.is_file():
                raise ValueError('Missing reference {}'.format(path))
            return str(path)
        check('FeGenie reference '+name, reference)
    if args.smoke_test and not failures:
        check('offline end-to-end smoke test', smoke_test)
    elif args.smoke_test:
        print('Smoke test skipped because prerequisites failed.', file=sys.stderr)
    if failures:
        print('Installation check failed: {}'.format(', '.join(failures)), file=sys.stderr)
        return 1
    print('Installation ready; no runtime downloads are required for bundled libraries.')
    return 0
