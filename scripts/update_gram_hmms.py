#!/usr/bin/env python3
"""Download the six official Pfam profiles; require --update for changed models."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from pathlib import Path
import sys
import tempfile
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from latest.models import profile_header

# Keep the interpretation attached to the versioned runtime manifest.
with (ROOT / "latest/manifests/GramGenie.tsv").open() as handle:
    RUNTIME = list(csv.DictReader(handle, delimiter="\t"))
MODELS = {r["model_id"]: (r["name"], r["description"], r["category"]) for r in RUNTIME}


def read_model(path, accession):
    h = profile_header(path)
    if h.get("ACC", "").split(".")[0] != accession:
        raise ValueError("Unexpected Pfam accession in {}".format(path))
    values = h.get("GA", "").rstrip(";").split()
    if len(values) != 2:
        raise ValueError("Missing sequence/domain GA cutoffs in {}".format(path))
    return {"versioned_accession": h["ACC"], "name": h["NAME"],
            "sequence_ga": float(values[0]), "domain_ga": float(values[1])}



def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=ROOT / "hmms" / "gram")
    p.add_argument("--update", action="store_true", help="Accept new model versions/checksums; review and commit changes")
    args = p.parse_args()
    destination = args.output
    existing = {}
    manifest = destination / "models.tsv"
    if manifest.exists():
        with manifest.open() as handle:
            existing = {r["accession"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with tempfile.TemporaryDirectory() as directory:
        staged, rows = Path(directory), []
        for acc, (name, description, role) in MODELS.items():
            url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{}/?annotation=hmm".format(acc)
            with urlopen(url, timeout=60) as response:
                content = response.read()
            if content.startswith(b"\x1f\x8b"):
                content = gzip.decompress(content)
            profile = staged / (acc + ".hmm")
            profile.write_bytes(content)
            metadata = read_model(profile, acc)
            digest = hashlib.sha256(content).hexdigest()
            if existing and not args.update and existing.get(acc, {}).get("sha256") != digest:
                raise SystemExit("{} changed upstream; use --update only after deciding to refresh the library".format(acc))
            rows.append({"accession": acc, "versioned_accession": metadata["versioned_accession"],
                         "name": metadata["name"], "description": description, "role": role,
                         "sequence_ga": metadata["sequence_ga"], "domain_ga": metadata["domain_ga"],
                         "sha256": digest, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "source_url": url})
        destination.mkdir(parents=True, exist_ok=True)
        for profile in staged.glob("*.hmm"):
            (destination / profile.name).write_bytes(profile.read_bytes())
        with manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    if destination.resolve() == (ROOT / "hmms/gram").resolve():
        for row in rows:
            runtime = next(r for r in RUNTIME if r["model_id"] == row["accession"])
            runtime.update(accession=row["versioned_accession"], name=row["name"],
                           sequence_cutoff=row["sequence_ga"], domain_cutoff=row["domain_ga"], sha256=row["sha256"])
        with (ROOT / "latest/manifests/GramGenie.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(RUNTIME[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(RUNTIME)
    print("Downloaded and validated all six Pfam profiles into {}".format(destination))


if __name__ == "__main__":
    main()
