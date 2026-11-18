"""Unit and optional real-HMMER integration checks (no network access required)."""
import csv
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
from genies_v3 import GramGenie as gram


class GramTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gram test ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_bundled_profiles_match_manifest(self):
        with (ROOT / "hmms/gram/models.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({r["accession"] for r in rows}, set(gram.MODELS))
        for row in rows:
            model = gram.read_model(ROOT / "hmms/gram" / (row["accession"] + ".hmm"), row["accession"])
            self.assertEqual(model["sha256"], row["sha256"])
            self.assertEqual(model["sequence_ga"], float(row["sequence_ga"]))
            self.assertEqual(model["domain_ga"], float(row["domain_ga"]))

    def test_fasta_identifiers_and_duplicate_rejection(self):
        path = self.root / "input.faa"
        path.write_text(">sp|P12345|name description\nMAAA*\n")
        self.assertEqual(gram.read_fasta(path), {"sp|P12345|name": "MAAA"})
        path.write_text(">same a\nMAAA\n>same b\nMCCC\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            gram.read_fasta(path)
        path.write_text("")
        with self.assertRaisesRegex(ValueError, "No protein"):
            gram.read_fasta(path)

    def test_sequence_and_domain_thresholds_and_repeat_retention(self):
        def row(score, domain, start):
            return "p - 200 LysM PF01476.28 44 1e-10 {} 0 1 3 1e-10 1e-10 {} 0 1 44 {} {} {} {} 0.99 desc\n".format(
                score, domain, start, start + 43, start, start + 43)
        path = self.root / "hits.domtblout"
        path.write_text(row(19, 12, 1) + row(30, 10, 1) + row(30, 12, 1) + row(30, 12, 80))
        models = {"PF01476": {"sequence_ga": 20.9, "domain_ga": 11.8, "name": "LysM"}}
        hits = gram.parse_domains(path, models, {"p": "M" * 200})
        self.assertEqual([h["ali_from"] for h in hits], [1, 80])

    def test_wrong_accession_or_missing_ga_rejected(self):
        source = (ROOT / "hmms/gram/PF01476.hmm").read_text()
        path = self.root / "bad.hmm"
        path.write_text(source)
        with self.assertRaisesRegex(ValueError, "Wrong"):
            gram.read_model(path, "PF04203")
        path.write_text("\n".join(line for line in source.splitlines() if not line.startswith("GA ")) + "\n")
        with self.assertRaisesRegex(ValueError, "GA cutoffs"):
            gram.read_model(path, "PF01476")

    def test_strict_dispatcher_arguments(self):
        r = subprocess.run([sys.executable, str(ROOT / "MagicLamp.py"), "GramGenie", "--not-a-real-option"],
                           capture_output=True, text=True, cwd=str(self.root))
        self.assertNotEqual(r.returncode, 0)

    def test_failed_subprocess_has_failed_manifest(self):
        inputs = self.root / "inputs"
        inputs.mkdir()
        (inputs / "sample.faa").write_text(">test\nMAAA\n")
        output = self.root / "output"
        with mock.patch.object(gram.shutil, "which", return_value="hmmsearch"), \
                mock.patch.object(gram.subprocess, "run", side_effect=subprocess.CalledProcessError(1, ["hmmsearch"])):
            self.assertEqual(gram.main(["-bin_dir", str(inputs), "-bin_ext", "faa", "--orfs", "-out", str(output)]), 1)
        self.assertEqual(json.loads((output / "gramgenie-run.json").read_text())["status"], "failed")
        self.assertFalse((output / "gramgenie.heatmap.csv").exists())

    @unittest.skipUnless(shutil.which("hmmsearch") and shutil.which("hmmemit"), "HMMER executables required")
    def test_real_hmmer_all_six_models_repeat_counts_and_zero_sample(self):
        inputs = self.root / "input with spaces"
        inputs.mkdir()
        sequences = {}
        for acc in gram.MODELS:
            result = subprocess.run(["hmmemit", "-c", str(ROOT / "hmms/gram" / (acc + ".hmm"))],
                                    check=True, capture_output=True, text=True)
            sequence = "".join(line for line in result.stdout.splitlines() if not line.startswith(">"))
            sequences[acc + "_control"] = (sequence + "GGGGSGGGGS") * 3
        # A multidomain protein must retain independent family annotations.
        sequences["multi_domain"] = sequences["PF01476_control"] + sequences["PF07501_control"]
        gram.write_fasta(inputs / "positive, sample.faa", sequences)
        (inputs / "negative.faa").write_text(">negative\n" + "A" * 150 + "\n")
        output = self.root / "output with spaces"
        env = dict(os.environ)
        env.pop("gram_hmms", None)
        r = subprocess.run([sys.executable, str(ROOT / "MagicLamp.py"), "GramGenie", "-bin_dir", str(inputs),
                            "-bin_ext", "faa", "--orfs", "--norm", "-out", str(output)],
                           cwd=str(self.root), env=env, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with (output / "gramgenie-summary.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        for acc in gram.MODELS:
            own = [r for r in rows if r["ORF"] == acc + "_control" and r["gene"] == acc]
            self.assertEqual(len(own), 1, acc)
            self.assertGreaterEqual(int(own[0]["domain_count"]), 2, acc)
        multi = {r["gene"] for r in rows if r["ORF"] == "multi_domain"}
        self.assertTrue({"PF01476", "PF07501"} <= multi)
        with (output / "gramgenie.heatmap.csv").open() as handle:
            heat = list(csv.reader(handle))
        zero_column = heat[0].index("negative.faa")
        self.assertTrue(all(float(row[zero_column]) == 0 for row in heat[1:]))
        self.assertEqual(len(heat), 7)
        self.assertEqual(len(rows), len({(r["file"], r["ORF"], r["gene"]) for r in rows}))
        denominator = len(sequences)
        positive_column = heat[0].index("positive, sample.faa")
        for row in heat[1:]:
            acc = row[0].split("_")[0]
            count = sum(r["gene"] == acc for r in rows)
            self.assertAlmostEqual(float(row[positive_column]), count * 1000 / denominator)
        self.assertEqual(json.loads((output / "gramgenie-run.json").read_text())["status"], "complete")
        # A rerun must not silently reuse old ORFs or overwrite existing results.
        self.assertEqual(gram.main(["-bin_dir", str(inputs), "-bin_ext", "faa", "--orfs", "-out", str(output)]), 1)

    @unittest.skipUnless(shutil.which("hmmsearch") and shutil.which("prodigal"), "HMMER and Prodigal required")
    def test_real_prodigal_zero_orfs_normalization(self):
        inputs = self.root / "nucleotides"
        inputs.mkdir()
        (inputs / "no_genes.fna").write_text(">contig\n" + "N" * 60 + "\n")
        output = self.root / "dna output"
        self.assertEqual(gram.main(["-bin_dir", str(inputs), "-bin_ext", "fna", "--meta", "--norm", "-out", str(output)]), 0)
        with (output / "gramgenie-samples.tsv").open() as handle:
            sample = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(sample["protein_count"], "0")
        with (output / "gramgenie.heatmap.csv").open() as handle:
            heat = list(csv.reader(handle))
        self.assertTrue(all(row[1] == "" for row in heat[1:]))


if __name__ == "__main__":
    unittest.main()
