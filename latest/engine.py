"""One checked execution path and output schema for every maintained genie."""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from . import __version__
from . import context, io
from .models import load_models, checksum, combined_profiles
from .registry import ROOT, GENIES, resolve_name


def make_parser(name):
    p = argparse.ArgumentParser(prog="MagicLamp.py " + name, description=GENIES[name][2])
    p.add_argument("-bin_dir", type=Path, required=True, help="Input directory")
    p.add_argument("-bin_ext", required=True, help="File extension, without the leading period")
    p.add_argument("-out", type=Path, default=Path(name.lower() + "_out"), help="New or empty output directory")
    p.add_argument("-t", type=int, default=1, help="HMMER threads per sample")
    group = p.add_mutually_exclusive_group()
    group.add_argument("--orfs", "--proteins", action="store_true", help="Protein FASTA; no genomic clustering")
    group.add_argument("--gbk", action="store_true", help="Annotated GenBank CDS translations; cluster by genomic order")
    p.add_argument("--meta", action="store_true", help="Prodigal metagenomic mode for nucleotide FASTA")
    p.add_argument("--no-cluster", action="store_true", help="Disable genomic neighborhoods for nucleotide/GenBank inputs")
    p.add_argument("-d", type=int, default=5, help="Maximum intervening CDSs between neighboring hits (default 5)")
    p.add_argument("-translation_table", type=int, default=11, help="Prodigal translation table")
    p.add_argument("--norm", action="store_true", help="Normalize counts per -inflation proteins")
    p.add_argument("-inflation", type=int, default=1000)
    p.add_argument("-hmm_dir", type=Path, help="Override the selected model-library root")
    p.add_argument("-cat", default="ALL", help="Restrict to a category (ALL by default)")
    p.add_argument("--makeplots", action="store_true", help="Make a count heatmap; requires matplotlib")
    p.add_argument("-ref", type=Path, help="Optional reference protein FASTA for separate BLASTP evidence")
    if name == "OmniGenie":
        p.add_argument("-genie", required=True, help="Supported library/name, e.g. iron, atp, gram")
    if name == "HmmGenie":
        p.add_argument("-eval", type=float, help="Explicit reporting E-value when custom profiles lack GA/TC")
        p.add_argument("--cutoff", choices=["auto", "ga", "tc", "evalue"], default="auto")
    if name in {"FeGenie", "OmniGenie"}:
        p.add_argument("--skip-reference", action="store_true", help="Explicitly skip bundled Geobacter/Thermincola BLASTP")
        p.add_argument("--heme", action="store_true", help="Report heme motifs across all proteins")
        p.add_argument("--hematite", action="store_true", help="Report legacy hematite-binding motif candidates")
    return p


def checked(command, log, manifest):
    manifest["commands"].append([str(x) for x in command])
    with Path(log).open("w") as handle:
        subprocess.run([str(x) for x in command], check=True, stdout=handle, stderr=subprocess.STDOUT)


def tool_version(executable):
    flag = "-version" if executable == "blastp" else "-v" if executable == "prodigal" else "-h"
    result = subprocess.run([executable, flag], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return " | ".join(lines[:3])


def search(models, combined, fasta, raw, tag, threads, lookup, manifest, classifier=False):
    model_by_key = {m.key: m for m in models}
    found = []
    for (policy, format_tag), path in combined.items():
        if (policy == "classifier") != classifier:
            continue
        search_tag = tag + "." + policy + "." + format_tag
        domfile = raw / (search_tag + ".domtblout")
        seqfile = raw / (search_tag + ".tblout")
        cmd = ["hmmsearch", "--cpu", threads, "--noali", "--domtblout", domfile, "--tblout", seqfile,
               "-o", raw / (search_tag + ".txt")]
        if policy in {"ga", "tc"}:
            cmd.append("--cut_" + policy)
        elif policy == "score":
            cmd += ["-E", .1, "--domE", 10]
        elif policy == "evalue":
            threshold = max(m.evalue for m in models if m.policy == policy)
            cmd += ["-E", threshold, "--domE", threshold]
        else:
            cmd += ["--max", "-T", -1000000, "--domT", -1000000]
        checked(cmd + [path, fasta], raw / (search_tag + ".log"), manifest)
        with domfile.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.split(maxsplit=22)
                if len(f) < 22 or f[3] not in model_by_key or f[0] not in lookup:
                    raise ValueError("Malformed/unrecognized HMMER domain row in {}".format(domfile))
                model = model_by_key[f[3]]
                full_e, full_score, domain_score = float(f[6]), float(f[7]), float(f[13])
                if model.sequence_cutoff is not None and full_score < model.sequence_cutoff:
                    continue
                if model.domain_cutoff is not None and domain_score < model.domain_cutoff:
                    continue
                if model.evalue is not None and full_e > model.evalue:
                    continue
                found.append({"target": f[0], "model": model, "evalue": full_e, "bit_score": full_score,
                              "domain_score": domain_score, "domain_evalue": float(f[12]),
                              "hmm_from": int(f[15]), "hmm_to": int(f[16]), "ali_from": int(f[17]),
                              "ali_to": int(f[18]), "env_from": int(f[19]), "env_to": int(f[20])})
    return found


def collapse(domains):
    grouped = defaultdict(list)
    for hit in domains:
        grouped[(hit["target"], hit["model"].model_id)].append(hit)
    hits = []
    for _, members in sorted(grouped.items()):
        best = dict(max(members, key=lambda h: h["domain_score"]))
        best["domain_count"] = len(members)
        hits.append(best)
    return hits


def ion_preferences(hits, models, paths, fasta, raw, tag, args, lookup, manifest):
    candidates = {h["target"] for h in hits if h["model"].accession.split(".")[0] in {"TIGR01260", "PF00137"}}
    if not candidates:
        return
    selected = {k: lookup[k] for k in sorted(candidates)}
    candidate_fasta = raw / (tag + ".proteolipids.faa")
    with candidate_fasta.open("w") as out:
        for k, protein in selected.items():
            out.write(">{}\n{}\n".format(k, protein.sequence))
    scores = defaultdict(dict)
    for hit in search(models, paths, candidate_fasta, raw, tag + ".ions", args.t, selected, manifest, classifier=True):
        scores[hit["target"]][hit["model"].name] = hit["bit_score"]
    for hit in hits:
        if hit["target"] not in candidates or hit["model"].accession.split(".")[0] not in {"TIGR01260", "PF00137"}:
            continue
        complex_type = hit["model"].properties.get("complex_type", "")
        prefix = "F-type" if complex_type == "F-type" else "V-type"
        local = {k: v for k, v in scores[hit["target"]].items() if k.startswith(prefix)}
        if not local:
            continue
        winner = max(local, key=local.get)
        others = [v for k, v in local.items() if k != winner]
        hit["ion_model"] = winner
        hit["ion_preference"] = "Na+" if "_Na_" in winner else "H+"
        # An unreported alternative is missing evidence, not a zero-bit score.
        hit["ion_score_delta"] = local[winner] - max(others) if others else ""
        hit["ion_confidence"] = "unresolved" if not others else (
            "high" if hit["ion_score_delta"] >= 20 else "medium" if hit["ion_score_delta"] >= 10 else
            "low" if hit["ion_score_delta"] >= 3 else "ambiguous")


def reference_hits(references, fasta, lookup, raw, tag, threads, manifest):
    results = []
    for index, reference in enumerate(references, 1):
        dest = raw / (tag + ".reference{}.tsv".format(index))
        cmd = ["blastp", "-query", fasta, "-subject", reference, "-evalue", "1e-10", "-num_threads", threads,
               "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore", "-out", dest]
        checked(cmd, raw / (tag + ".reference{}.log".format(index)), manifest)
        with dest.open() as handle:
            for line in handle:
                f = line.rstrip().split("\t")
                if len(f) != 10 or f[0] not in lookup:
                    raise ValueError("Unrecognized BLAST row in {}".format(dest))
                results.append({"target": f[0], "reference_library": reference.name, "reference_id": f[1],
                                "percent_identity": float(f[2]), "alignment_length": int(f[3]),
                                "query_start": int(f[4]), "query_end": int(f[5]), "reference_start": int(f[6]),
                                "reference_end": int(f[7]), "evalue": float(f[8]), "bit_score": float(f[9])})
    return results


COMMON = ["sample", "protein_id", "model", "accession", "hmm_name", "gene", "description", "category",
          "annotation_source", "evalue", "bit_score", "cutoff_policy", "sequence_cutoff", "domain_cutoff",
          "domain_count", "contig", "start", "end", "strand", "gene_order", "cluster_id", "context_status",
          "context_annotation", "complex_completeness", "ion_model", "ion_preference", "ion_score_delta",
          "ion_confidence", "heme_c_motifs", "heme_b_motifs", "sequence"]
DOMAIN = ["sample", "protein_id", "model", "accession", "domain_score", "domain_evalue", "hmm_from", "hmm_to",
          "ali_from", "ali_to", "env_from", "env_to"]


def table(path, fields, rows, delimiter="\t"):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(out, name, models, samples, rows, domains, clusters, references, args, motifs):
    table(out / "summary.csv", COMMON, rows, ",")
    table(out / "domains.tsv", DOMAIN, domains)
    table(out / "samples.tsv", ["sample", "protein_count", "proteins_with_hit", "domain_hits", "status"], samples)
    table(out / "clusters.tsv", ["sample", "cluster_id", "contig", "start", "end", "protein_count", "proteins"], clusters)
    table(out / "reference_hits.tsv", ["sample", "protein_id", "reference_library", "reference_id", "percent_identity",
                                       "alignment_length", "query_start", "query_end", "reference_start", "reference_end",
                                       "evalue", "bit_score"], references)
    if getattr(args, "heme", False) or getattr(args, "hematite", False):
        table(out / "motif_hits.tsv", ["sample", "protein_id", "heme_c_motifs", "heme_b_motifs", "hematite_motifs", "sequence"], motifs)
    primary_models = [m for m in models if m.policy != "classifier"]
    labels = [m.model_id for m in primary_models]
    categories = sorted({m.category for m in primary_models})
    model_counts, category_counts = defaultdict(set), defaultdict(set)
    for row in rows:
        model_counts[(row["sample"], row["model"])].add(row["protein_id"])
        category_counts[(row["sample"], row["category"])].add(row["protein_id"])
    for filename, values, count_dict in [("counts.csv", labels, model_counts), ("categories.csv", categories, category_counts)]:
        with (out / filename).open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["X"] + [s["sample"] for s in samples])
            for label in values:
                counts = []
                for sample in samples:
                    count = len(count_dict[(sample["sample"], label)])
                    counts.append((count * args.inflation / sample["protein_count"] if sample["protein_count"] else "")
                                  if args.norm else count)
                writer.writerow([label] + counts)
    if name == "FeGenie":
        supported = [r for r in rows if r["context_status"] in {"supported_neighborhood", "sequence_supported_candidate"}]
        table(out / "context_supported.csv", COMMON, supported, ",")
    # Stable aliases keep the familiar per-genie naming convention.
    shutil.copyfile(out / "summary.csv", out / (name.lower() + "-summary.csv"))
    shutil.copyfile(out / "counts.csv", out / (name.lower() + ".heatmap.csv"))
    if args.makeplots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        matrix = np.array([[len(model_counts[(s["sample"], label)]) for s in samples] for label in labels], dtype=float)
        if args.norm:
            denom = np.array([s["protein_count"] for s in samples], dtype=float)
            matrix = np.divide(matrix * args.inflation, denom, out=np.full_like(matrix, np.nan), where=denom > 0)
        # Plotting is a view of the same matrix, never a separate normalization.
        fig, ax = plt.subplots(figsize=(max(6, min(20, len(samples))), max(3, min(30, len(labels) * .18))))
        plot = ax.imshow(matrix, aspect="auto", interpolation="none")
        ax.set_xticks(range(len(samples)), [s["sample"] for s in samples], rotation=90)
        if len(labels) <= 100:
            ax.set_yticks(range(len(labels)), labels)
        ax.set_xlabel("Sample")
        ax.set_ylabel("HMM model")
        fig.colorbar(plot, ax=ax, label="Proteins per {}".format(args.inflation) if args.norm else "Protein count")
        fig.tight_layout()
        fig.savefig(out / "heatmap.png", dpi=150)
        plt.close(fig)


def run(name, args):
    if args.t < 1 or args.d < 0 or args.inflation < 1:
        raise ValueError("Threads/inflation must be positive and -d must be nonnegative")
    if args.meta and (args.orfs or args.gbk):
        raise ValueError("--meta applies only to nucleotide FASTA")
    if not args.bin_dir.is_dir():
        raise ValueError("Input directory does not exist: {}".format(args.bin_dir))
    extension = args.bin_ext.lstrip(".")
    if not extension or "/" in extension or "\\" in extension:
        raise ValueError("Invalid input file extension")
    inputs = sorted(p for p in args.bin_dir.iterdir() if p.is_file() and p.name.endswith("." + extension))
    if not inputs:
        raise ValueError("No matching input files")
    if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
        raise ValueError("Output must be a new or empty directory: {}".format(args.out))
    effective = resolve_name(args.genie) if name == "OmniGenie" else name
    if effective in {"HmmGenie", "OmniGenie"} and name == "OmniGenie":
        raise ValueError("OmniGenie requires a named supported library")
    evalue = getattr(args, "eval", None)
    if evalue is not None and (not math.isfinite(evalue) or evalue <= 0):
        raise ValueError("-eval must be a finite positive value")
    models, library = load_models(effective, args.hmm_dir, evalue, getattr(args, "cutoff", "auto"))
    if args.cat != "ALL":
        selected = [m for m in models if m.policy != "classifier" and
                    (m.category == args.cat or m.properties.get("reaction") == args.cat)]
        if not selected:
            raise ValueError("Unknown/empty category {!r}".format(args.cat))
        models = selected + [m for m in models if m.policy == "classifier"]
    references = [args.ref.expanduser().resolve()] if args.ref else []
    if effective == "FeGenie" and not getattr(args, "skip_reference", False):
        references += [library / "iron_reduction/non-aligned" / f for f in ["TherJR_SLCs.faa", "geobacter_PCCs.faa"]]
    for reference in references:
        if not reference.is_file():
            raise ValueError("Reference FASTA missing: {}. FeGenie can explicitly use --skip-reference.".format(reference))
    needed = ["hmmsearch"] + ([] if args.orfs or args.gbk else ["prodigal"]) + (["blastp"] if references else [])
    for executable in needed:
        if not shutil.which(executable):
            raise ValueError("Required executable missing from PATH: {}".format(executable))
    if args.gbk:
        try:
            import Bio
        except ImportError as exc:
            raise ValueError("--gbk requires Biopython") from exc
    if args.makeplots:
        try:
            import matplotlib, numpy
        except ImportError as exc:
            raise ValueError("--makeplots requires matplotlib and numpy (see environment-plots.yml)") from exc
    cluster_enabled = not args.orfs and not args.no_cluster
    args.out.mkdir(parents=True, exist_ok=True)
    out = args.out.resolve()
    raw, orfs = out / "HMM_results", out / "ORF_calls"
    raw.mkdir()
    orfs.mkdir()
    manifest = {"suite_version": __version__, "genie": name, "effective_library": effective, "status": "running",
                "started_utc": datetime.now(timezone.utc).isoformat(), "input_format": "orfs" if args.orfs else "gbk" if args.gbk else "nucleotide",
                "clustering": cluster_enabled, "max_intervening_cds": args.d if cluster_enabled else None,
                "normalized_counts": args.norm, "inflation": args.inflation,
                "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
                "python": sys.version, "tools": {}, "python_packages": {},
                "models": [m.provenance() for m in models], "samples": [], "commands": [],
                "references": [{"path": str(p), "sha256": checksum(p)} for p in references]}
    def save_manifest():
        payload = json.dumps(manifest, indent=2) + "\n"
        (out / "run.json").write_text(payload)
        (out / (name.lower() + "-run.json")).write_text(payload)
    save_manifest()
    samples, rows, all_domains, all_clusters, all_references, motifs = [], [], [], [], [], []
    try:
        manifest["tools"] = {exe: tool_version(exe) for exe in needed}
        packages = (["biopython"] if args.gbk else []) + (["matplotlib", "numpy"] if args.makeplots else [])
        manifest["python_packages"] = {package: importlib.metadata.version(package) for package in packages}
        paths = combined_profiles(models, raw)
        for index, source in enumerate(inputs, 1):
            tag = "sample{:06d}".format(index)
            input_hash = checksum(source)
            if args.orfs:
                proteins = io.protein_fasta(source)
            elif args.gbk:
                proteins = io.genbank(source)
            else:
                predicted = orfs / (tag + ".prodigal.faa")
                cmd = ["prodigal", "-i", source.resolve(), "-a", predicted, "-o", orfs / (tag + ".prodigal.gbk"),
                       "-g", args.translation_table, "-p", "meta" if args.meta else "single", "-q"]
                checked(cmd, orfs / (tag + ".prodigal.log"), manifest)
                proteins = io.prodigal_fasta(predicted)
            fasta = orfs / (tag + ".faa")
            lookup = io.write_search_fasta(fasta, proteins)
            table(orfs / (tag + ".protein_ids.tsv"), ["search_id", "protein_id", "contig", "start", "end", "strand", "gene_order"],
                  [{"search_id": k, "protein_id": p.identifier, "contig": p.contig, "start": p.start, "end": p.end,
                    "strand": p.strand, "gene_order": p.order} for k, p in lookup.items()])
            domains = search(models, paths, fasta, raw, tag, args.t, lookup, manifest) if lookup else []
            hits = collapse(domains)
            clusters = context.assign_clusters(hits, lookup, cluster_enabled, args.d)
            context.annotate(effective, hits, lookup, cluster_enabled)
            if effective == "ATPGenie" and lookup:
                ion_preferences(hits, models, paths, fasta, raw, tag, args, lookup, manifest)
            refs = reference_hits(references, fasta, lookup, raw, tag, args.t, manifest) if lookup and references else []
            for ref in refs:
                ref["sample"], ref["protein_id"] = source.name, lookup[ref.pop("target")].identifier
                all_references.append(ref)
            for hit in hits:
                model, protein = hit["model"], lookup[hit["target"]]
                row = {key: hit.get(key, "") for key in COMMON}
                row.update(sample=source.name, protein_id=protein.identifier, model=model.model_id,
                           accession=model.accession, hmm_name=model.name, gene=model.gene, description=model.description,
                           category=model.category, annotation_source=model.annotation_source, cutoff_policy=model.policy,
                           sequence_cutoff=model.sequence_cutoff if model.sequence_cutoff is not None else "",
                           domain_cutoff=model.domain_cutoff if model.domain_cutoff is not None else "",
                           contig=protein.contig, start=protein.start, end=protein.end, strand=protein.strand,
                           gene_order=protein.order, sequence=protein.sequence)
                rows.append(row)
            for hit in domains:
                all_domains.append(dict(sample=source.name, protein_id=lookup[hit["target"]].identifier,
                                        model=hit["model"].model_id, accession=hit["model"].accession,
                                        **{k: hit[k] for k in DOMAIN[4:] if k in hit}))
            for cluster in clusters:
                all_clusters.append(dict(sample=source.name, **cluster))
            sample = {"sample": source.name, "protein_count": len(lookup), "proteins_with_hit": len({h["target"] for h in hits}),
                      "domain_hits": len(domains), "status": "searched" if lookup else "no_translated_proteins"}
            samples.append(sample)
            manifest["samples"].append(dict(path=str(source.resolve()), sha256=input_hash, search_fasta=str(fasta), **sample))
            if getattr(args, "heme", False) or getattr(args, "hematite", False):
                for protein in proteins:
                    hc = context.heme_motifs(protein.sequence)
                    hb = len(re.findall(r"G.[HR].C[PLAV]G", protein.sequence))
                    hm = len(re.findall(r"[STC][AVILMFYWH][ST]P[ST]", protein.sequence)) if getattr(args, "hematite", False) else ""
                    if (getattr(args, "heme", False) and (hc or hb)) or hm:
                        motifs.append({"sample": source.name, "protein_id": protein.identifier, "heme_c_motifs": hc,
                                       "heme_b_motifs": hb, "hematite_motifs": hm, "sequence": protein.sequence})
        write_outputs(out, name, models, samples, rows, all_domains, all_clusters, all_references, args, motifs)
        manifest["status"] = "complete"
    except BaseException as exc:
        manifest["status"], manifest["error"] = "failed", str(exc)
        raise
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save_manifest()
    print("{} results: {}".format(name, out))


def main(name, argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == name:
        arguments.pop(0)
    args = make_parser(name).parse_args(arguments)
    try:
        run(name, args)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print("{}: {}".format(name, exc), file=sys.stderr)
        return 1
    return 0
