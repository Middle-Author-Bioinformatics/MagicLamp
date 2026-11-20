"""Installation metadata and dependency diagnostics, independent of network access."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tomllib
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from latest import __version__
from latest import doctor


class InstallationTests(unittest.TestCase):
    def test_package_version_and_entry_points(self):
        with (ROOT/'pyproject.toml').open('rb') as handle:
            config = tomllib.load(handle)
        self.assertEqual(config['project']['version'], __version__)
        self.assertEqual(config['project']['scripts']['magiclamp'], 'latest.cli:main')
        self.assertEqual(config['project']['scripts']['MagicLamp.py'], 'latest.cli:main')
        self.assertEqual(set(config['tool']['setuptools']['packages']), {'latest','magiclamp_models'})

    def test_fixture_matches_previous_smoke_test(self):
        evidence = json.loads((ROOT/'docs/audit/latest_nucleotide_smoke.json').read_text())
        digest = hashlib.sha256((ROOT/'latest/fixtures/B_subtilis_first_30kb.fna').read_bytes()).hexdigest()
        self.assertEqual(digest, evidence['input']['sha256'])

    def test_external_tool_version_detection(self):
        for tool, output, minimum in [('hmmsearch','# HMMER 3.4 (Aug 2023)',(3,3)),
                                      ('prodigal','Prodigal V2.6.3: February, 2016',(2,6,3)),
                                      ('blastp','blastp: 2.17.0+',(2,12))]:
            with self.subTest(tool=tool), mock.patch.object(doctor.shutil,'which',return_value='/tools/'+tool), \
                 mock.patch.object(doctor.subprocess,'run',return_value=subprocess.CompletedProcess([],0,stdout=output)):
                self.assertIn('/tools/'+tool, doctor.check_tool(tool, [], minimum))

    def test_old_or_missing_tool_fails(self):
        with mock.patch.object(doctor.shutil,'which',return_value=None):
            with self.assertRaisesRegex(ValueError,'missing from PATH'):
                doctor.check_tool('hmmsearch', ['-h'], (3,3))
        with mock.patch.object(doctor.shutil,'which',return_value='/tools/hmmsearch'), \
             mock.patch.object(doctor.subprocess,'run',return_value=subprocess.CompletedProcess([],0,stdout='HMMER 3.1b2')):
            with self.assertRaisesRegex(ValueError,'too old'):
                doctor.check_tool('hmmsearch', ['-h'], (3,3))

    def test_doctor_does_not_run_controls_after_failed_prerequisite(self):
        with mock.patch.object(doctor,'check_tool',side_effect=ValueError('missing executable')), \
             mock.patch.object(doctor,'smoke_test') as smoke, \
             mock.patch('sys.stdout'), mock.patch('sys.stderr'):
            self.assertEqual(doctor.main(['--smoke-test']),1)
            smoke.assert_not_called()


if __name__ == '__main__':
    unittest.main()
