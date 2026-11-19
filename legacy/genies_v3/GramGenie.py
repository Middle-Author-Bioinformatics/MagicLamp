#!/usr/bin/env python3
"""Profile six Pfam cell-wall-associated domains with their gathering cutoffs."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


MODELS = {
    "PF01473": ("Choline_bind_1", "Putative cell wall binding repeat", "cell_wall_binding"),
    "PF04122": ("CW_binding_2", "Cell wall binding domain 2 (CWB2)", "cell_wall_binding"),
    "PF00395": ("SLH", "S-layer homology domain", "cell_wall_binding"),
    "PF04203": ("Sortase", "Sortase domain", "anchoring_enzyme"),
    "PF01476": ("LysM", "LysM domain", "cell_wall_associated"),
    "PF07501": ("G5", "G5 domain", "cell_wall_associated"),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def read_model(path, expected):
    text = Path(path).read_text()
    if not text.startswith("HMMER3/") or text.count("\nNAME ") != 1 or text.count("\n//") != 1:
        raise ValueError("Expected exactly one HMMER3 profile in {}".format(path))
    fields = {}
    for line in text.splitlines():
        if line.startswith("HMM "):
            break
        parts = line.split(None, 1)
        if len(parts) == 2:
            fields[parts[0]] = parts[1]
    if fields.get("ACC", "").split(".")[0] != expected:
        raise ValueError("Wrong or missing ACC in {} (expected {})".format(path, expected))
    if fields.get("ALPH") != "amino":
        raise ValueError("Expected an amino-acid HMM: {}".format(path))
    try:
        sequence_ga, domain_ga = map(float, fields["GA"].rstrip(";").split())
    except (KeyError, ValueError):
        raise ValueError("Missing/invalid Pfam sequence and domain GA cutoffs: {}".format(path))
    return {
        "accession": expected, "versioned_accession": fields["ACC"],
        "name": fields["NAME"], "length": int(fields["LENG"]),
        "sequence_ga": sequence_ga, "domain_ga": domain_ga,
        "sha256": sha256(path), "path": str(Path(path).resolve()),
    }


def read_fasta(path, allow_empty=False):
    records = {}
    header, sequence = None, []

    def finish():
        if header is None:
            return
        seq = "".join(sequence).upper().rstrip("*")
        if not seq or set(seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError("Empty or invalid protein sequence {} in {}".format(header, path))
        if header in records:
            raise ValueError("Duplicate protein ID {} in {}".format(header, path))
        records[header] = seq

    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                parts = line[1:].split()
                if not parts:
                    raise ValueError("Empty FASTA identifier in {}".format(path))
                header, sequence = parts[0], []
            else:
                if header is None:
                    raise ValueError("Sequence before FASTA header in {}".format(path))
                sequence.append(line)
    finish()
    if not records and not allow_empty:
        raise ValueError("No protein sequences found in {}".format(path))
    return records


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for name, sequence in records.items():
            handle.write(">{}\n{}\n".format(name, sequence))


def genbank_proteins(path):
    try:
        from Bio import SeqIO
    except ImportError:
        raise ValueError("--gbk requires Biopython; alternatively supply translated FASTA with --orfs")
    records = {}
    for record in SeqIO.parse(str(path), "genbank"):
        for index, feature in enumerate(record.features, 1):
            if feature.type != "CDS" or "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
                continue
            translation = feature.qualifiers.get("translation", [])
            if not translation:
                continue
            tag = feature.qualifiers.get("locus_tag", feature.qualifiers.get("protein_id", [str(index)]))[0]
            name = "{}|{}".format(record.id, tag)
            if name in records:
                raise ValueError("Duplicate CDS identifier {} in {}".format(name, path))
            records[name] = translation[0]
    if not records:
        raise ValueError("No translated, non-pseudogene CDS features in {}".format(path))
    return records


def parse_domains(path, models, proteins):
    hits = []
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError("Malformed domtblout row {} in {}".format(number, path))
            accession = fields[4].split(".")[0]
            if accession not in models or fields[0] not in proteins:
                raise ValueError("Unexpected model or protein in {}:{}".format(path, number))
            model = models[accession]
            # Check both thresholds explicitly, in addition to hmmsearch --cut_ga.
            if float(fields[7]) < model["sequence_ga"] or float(fields[13]) < model["domain_ga"]:
                continue
            hits.append({
                "ORF": fields[0], "gene": accession, "hmm_name": model["name"],
                "hmm_accession": fields[4], "protein_length": int(fields[2]),
                "hmm_length": int(fields[5]), "sequence_evalue": float(fields[6]),
                "sequence_score": float(fields[7]), "domain_number": int(fields[9]),
                "independent_evalue": float(fields[12]), "domain_score": float(fields[13]),
                "hmm_from": int(fields[15]), "hmm_to": int(fields[16]),
                "ali_from": int(fields[17]), "ali_to": int(fields[18]),
                "env_from": int(fields[19]), "env_to": int(fields[20]),
            })
    return sorted(hits, key=lambda h: (h["ORF"], h["gene"], h["ali_from"]))


def parser():
    p = argparse.ArgumentParser(
        prog="MagicLamp.py GramGenie",
        description="Detect six Pfam cell-wall-associated domains using sequence and domain GA cutoffs. "
                    "Hits do not establish Gram-positive identity, localization, or EET activity.",
    )
    p.add_argument("-bin_dir", required=True, type=Path, help="Input directory")
    p.add_argument("-bin_ext", required=True, help="Input extension (e.g. fna or faa)")
    p.add_argument("-out", type=Path, default=Path("gramgenie_out"), help="New or empty output directory")
    p.add_argument("-t", type=int, default=1, help="HMMER CPU threads per sample")
    formats = p.add_mutually_exclusive_group()
    formats.add_argument("--orfs", "--proteins", action="store_true", help="Input is protein FASTA; skip Prodigal")
    formats.add_argument("--gbk", action="store_true", help="Extract annotated CDS translations (requires Biopython)")
    p.add_argument("--meta", action="store_true", help="Use Prodigal metagenome mode for nucleotide FASTA")
    p.add_argument("--norm", action="store_true", help="Normalize heatmap counts per -inflation proteins")
    p.add_argument("-inflation", type=int, default=1000)
    p.add_argument("-translation_table", type=int, default=11, help="Prodigal translation table")
    p.add_argument("-hmm_dir", type=Path, help="Override $gram_hmms or bundled hmms/gram")
    return p


def run(args):
    if args.t < 1 or args.inflation < 1:
        raise ValueError("-t and -inflation must be positive integers")
    if args.meta and (args.orfs or args.gbk):
        raise ValueError("--meta applies only to nucleotide FASTA")
    if not args.bin_dir.is_dir():
        raise ValueError("Input directory does not exist: {}".format(args.bin_dir))
    extension = args.bin_ext.lstrip(".")
    if not extension or "/" in extension or "\\" in extension:
        raise ValueError("Provide a file extension such as fna or faa")
    inputs = sorted(p for p in args.bin_dir.iterdir() if p.is_file() and p.name.endswith("." + extension))
    if not inputs:
        raise ValueError("No .{} input files in {}".format(extension, args.bin_dir))
    if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
        raise ValueError("Output must be a new or empty directory: {}".format(args.out))
    library = args.hmm_dir or Path(os.environ.get("gram_hmms") or Path(__file__).resolve().parents[1] / "hmms" / "gram")
    models = {acc: read_model(library / (acc + ".hmm"), acc) for acc in MODELS}
    if not shutil.which("hmmsearch"):
        raise ValueError("hmmsearch is required (install HMMER and put it on PATH)")
    if not args.orfs and not args.gbk and not shutil.which("prodigal"):
        raise ValueError("prodigal is required for nucleotide FASTA; use --orfs for proteins")

    args.out.mkdir(parents=True, exist_ok=True)
    orfs, raw = args.out / "ORF_calls", args.out / "HMM_results"
    orfs.mkdir()
    raw.mkdir()
    manifest = {
        "genie": "GramGenie", "version": 1, "status": "running", "cutoff": "Pfam GA (sequence and domain)",
        "models": list(models.values()), "normalized_heatmap": args.norm, "inflation": args.inflation,
        "input_format": "proteins" if args.orfs else ("genbank" if args.gbk else "nucleotide"),
        "translation_table": args.translation_table if not (args.orfs or args.gbk) else None,
        "prodigal_mode": "meta" if args.meta else "single", "commands": [], "samples": [],
    }
    manifest_path = args.out / "gramgenie-run.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    combined = raw / "GramGenie.hmm"
    combined.write_text("".join((library / (acc + ".hmm")).read_text().rstrip() + "\n" for acc in MODELS))
    all_domains, summary, sample_info, counts = [], [], [], {}
    try:
        for source in inputs:
            proteins_path = orfs / (source.name + "-proteins.faa")
            if args.orfs:
                proteins = read_fasta(source)
            elif args.gbk:
                write_fasta(proteins_path, genbank_proteins(source))
                proteins = read_fasta(proteins_path)
            else:
                command = ["prodigal", "-i", str(source.resolve()), "-a", str(proteins_path),
                           "-o", str(orfs / (source.name + "-prodigal.gbk")), "-g", str(args.translation_table),
                           "-p", "meta" if args.meta else "single", "-q"]
                manifest["commands"].append(command)
                with (orfs / (source.name + "-prodigal.log")).open("w") as log:
                    subprocess.run(command, check=True, stdout=log, stderr=log)
                proteins = read_fasta(proteins_path, allow_empty=True)
            # Normalize identifiers and terminal stop symbols consistently before HMMER.
            write_fasta(proteins_path, proteins)
            sample = {"file": source.name, "path": str(source.resolve()), "sha256": sha256(source),
                      "protein_count": len(proteins)}
            manifest["samples"].append(sample)
            domain_path = raw / (source.name + ".domtblout")
            if proteins:
                command = ["hmmsearch", "--cpu", str(args.t), "--cut_ga", "--noali", "--domtblout", str(domain_path),
                           "-o", str(raw / (source.name + ".hmmsearch.txt")), str(combined), str(proteins_path)]
                manifest["commands"].append(command)
                with (raw / (source.name + ".hmmsearch.log")).open("w") as log:
                    subprocess.run(command, check=True, stdout=log, stderr=log)
                domains = parse_domains(domain_path, models, proteins)
            else:
                domain_path.write_text("# No predicted proteins; no search performed.\n")
                domains = []
            counts[source.name] = {acc: set() for acc in MODELS}
            grouped = {}
            for hit in domains:
                all_domains.append(dict(file=source.name, **hit))
                counts[source.name][hit["gene"]].add(hit["ORF"])
                grouped.setdefault((hit["ORF"], hit["gene"]), []).append(hit)
            for (protein, accession), hits in sorted(grouped.items()):
                name, description, category = MODELS[accession]
                best = max(hits, key=lambda h: h["domain_score"])
                summary.append({"file": source.name, "ORF": protein, "gene": accession, "hmm_name": name,
                                "description": description, "subcategory": category,
                                "evalue": best["sequence_evalue"], "bit_score": best["sequence_score"],
                                "domain_count": len(hits), "seq": proteins[protein]})
            sample_info.append({"file": source.name, "protein_count": len(proteins),
                                "proteins_with_hit": len({h["ORF"] for h in domains}),
                                "domain_hits": len(domains)})

        def table(filename, fields, rows, delimiter=","):
            with (args.out / filename).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
                writer.writeheader()
                writer.writerows(rows)

        table("gramgenie-summary.csv", ["file", "ORF", "gene", "hmm_name", "description", "subcategory",
                                         "evalue", "bit_score", "domain_count", "seq"], summary)
        table("gramgenie-domains.tsv", ["file", "ORF", "gene", "hmm_name", "hmm_accession", "protein_length",
                                        "hmm_length", "sequence_evalue", "sequence_score", "domain_number",
                                        "independent_evalue", "domain_score", "hmm_from", "hmm_to", "ali_from",
                                        "ali_to", "env_from", "env_to"], all_domains, "\t")
        table("gramgenie-samples.tsv", ["file", "protein_count", "proteins_with_hit", "domain_hits"], sample_info, "\t")
        with (args.out / "gramgenie.heatmap.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["X"] + [p.name for p in inputs])
            denominators = {r["file"]: r["protein_count"] for r in sample_info}
            for acc, meta in MODELS.items():
                values = []
                for source in inputs:
                    count = len(counts[source.name][acc])
                    denominator = denominators[source.name]
                    values.append((count * args.inflation / denominator if denominator else "") if args.norm else count)
                writer.writerow([acc + "_" + meta[0]] + values)
        manifest["status"] = "complete"
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print("GramGenie results: {}".format(args.out))


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "GramGenie":
        arguments.pop(0)
    p = parser()
    args = p.parse_args(arguments)
    try:
        run(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print("GramGenie: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
