"""Validated, versioned model manifests and explicit threshold policies."""
import csv
from dataclasses import dataclass, field
import hashlib
import json
import math
import os
from pathlib import Path
import re
from .registry import HMM_ROOT, GENIES


def checksum(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def profile_header(path):
    text = Path(path).read_text()
    if not text.startswith("HMMER3/") or len(re.findall(r"^NAME\s", text, re.M)) != 1 or text.count("\n//") != 1:
        raise ValueError("Expected one complete HMMER3 profile: {}".format(path))
    fields = {}
    for line in text.splitlines():
        if line.startswith("HMM "):
            break
        parts = line.split(None, 1)
        if len(parts) == 2:
            fields[parts[0]] = parts[1]
    if fields.get("ALPH") != "amino":
        raise ValueError("Expected an amino-acid HMM: {}".format(path))
    return fields


@dataclass
class Model:
    model_id: str
    path: Path
    accession: str
    name: str
    gene: str
    description: str
    category: str
    policy: str
    sequence_cutoff: object
    domain_cutoff: object
    evalue: object
    digest: str
    annotation_source: str
    properties: dict = field(default_factory=dict)
    key: str = ""

    def provenance(self):
        return {"model": self.model_id, "path": str(self.path), "accession": self.accession,
                "name": self.name, "gene": self.gene, "description": self.description,
                "category": self.category, "policy": self.policy, "properties": self.properties,
                "sequence_cutoff": self.sequence_cutoff, "domain_cutoff": self.domain_cutoff,
                "evalue": self.evalue, "sha256": self.digest, "annotation_source": self.annotation_source}


def number(value):
    if value in ("", None):
        return None
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Non-finite model cutoff")
    return value


def load_models(name, override=None, evalue=None, cutoff="auto"):
    if name == "HmmGenie":
        if override is None:
            raise ValueError("HmmGenie requires -hmm_dir")
        library = Path(override).expanduser().resolve()
        paths = sorted(p for p in library.rglob("*") if p.is_file() and p.suffix.lower() == ".hmm")
        models = []
        for path in paths:
            h = profile_header(path)
            policy = cutoff if cutoff != "auto" else ("ga" if "GA" in h else "tc" if "TC" in h else "evalue")
            seq, dom = None, None
            if policy in {"ga", "tc"}:
                values = h.get(policy.upper(), "").rstrip(";").split()
                if len(values) != 2:
                    raise ValueError("No {} cutoffs in {}".format(policy.upper(), path))
                seq, dom = map(number, values)
            elif evalue is None:
                raise ValueError("No GA/TC thresholds in {}; explicitly supply -eval for this custom library".format(path))
            models.append(Model(str(path.relative_to(library).with_suffix("")), path, h.get("ACC", ""), h["NAME"],
                                h["NAME"], h.get("DESC", h["NAME"]), "custom", policy, seq, dom, evalue,
                                checksum(path), "custom HMM header"))
    else:
        directory, variable, _ = GENIES[name]
        library = Path(override or os.environ.get(variable) or HMM_ROOT / directory).expanduser().resolve()
        manifest = Path(__file__).parent / "manifests" / (name + ".tsv")
        with manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        models = []
        for row in rows:
            path = (library / row["path"]).resolve()
            if not path.is_relative_to(library):
                raise ValueError("Model path escapes its library: {}".format(path))
            h = profile_header(path)
            digest = checksum(path)
            if digest != row["sha256"]:
                raise ValueError("Model checksum changed: {}. Review/update the versioned manifest.".format(path))
            policy = row["cutoff_policy"]
            seq, dom = number(row["sequence_cutoff"]), number(row["domain_cutoff"])
            if policy in {"tc", "ga"}:
                actual = [number(v) for v in h.get(policy.upper(), "").rstrip(";").split()]
                if actual != [seq, dom]:
                    raise ValueError("Manifest cutoff does not match HMM: {}".format(path))
            if policy == "score" and seq is None:
                raise ValueError("Missing curated cutoff: {}".format(path))
            models.append(Model(row["model_id"], path, row["accession"], row["name"], row["gene"],
                                row["description"], row["category"], policy, seq, dom, number(row["evalue"]),
                                digest, row["annotation_source"], json.loads(row["properties"])))
    if not models:
        raise ValueError("No validated HMM models in {}".format(library))
    if len({m.model_id for m in models}) != len(models):
        raise ValueError("Duplicate model identifiers in {}".format(library))
    for i, model in enumerate(models, 1):
        model.key = "m{:06d}".format(i)
    return models, library


def combined_profiles(models, directory):
    groups = {}
    for model in models:
        # HMMER selects a parser from the first record. Mixing 3/b and 3/f
        # in one file fails even though both profiles are valid individually.
        with model.path.open() as handle:
            format_tag = handle.readline().split()[0].replace("/", "")
        groups.setdefault((model.policy, format_tag), []).append(model)
    paths = {}
    for (policy, format_tag), members in groups.items():
        path = directory / (policy + "_" + format_tag + ".hmm")
        with path.open("w") as handle:
            for model in members:
                # Query IDs must be unambiguous even when source NAME/ACC fields repeat.
                text = re.sub(r"^NAME\s+.*$", "NAME  " + model.key, model.path.read_text(), count=1, flags=re.M)
                handle.write(text.rstrip() + "\n")
        paths[(policy, format_tag)] = path
    return paths
