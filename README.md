# MagicLamp — latest

A shared, maintained suite for screening protein sequences and genomes against
curated HMM libraries. `MagicLamp.py` dispatches only the `latest/` package
(version 4.0.1). Historical implementations are labeled and retained in `legacy/`.
AbxGenie and BTEXGenie are excluded from this version.

## Install

[Complete installation and GramGenie walkthrough](docs/INSTALL.md) includes
Miniforge setup, prerequisites, input formats, outputs, upgrades, and an offline
installation test.

```bash
git clone https://github.com/Middle-Author-Bioinformatics/MagicLamp.git
cd MagicLamp
bash setup.sh --name magiclamp --plots
conda activate magiclamp
magiclamp doctor --smoke-test --plots

# One assembled nucleotide FASTA per genome; use your own paths.
magiclamp GramGenie -bin_dir /path/to/genomes -bin_ext fna \
  -out /path/to/gramgenie_results -t 24 --makeplots
```

The installer creates a new environment and installs a normal package including
its HMM profiles/manifests and FeGenie reference sequences. It does not require
keeping the clone in place. The commands `magiclamp` and `MagicLamp.py` both work
from any directory; direct `python MagicLamp.py` execution in a checkout remains
supported. Bundled analyses work offline after setup.

The core environment declares Python 3.11, HMMER, Prodigal, BLAST+, Biopython,
pip, setuptools, and wheel. `--plots` adds matplotlib/numpy. HMMER, Prodigal, and
BLAST+ are external programs; pip alone cannot install them. R, DIAMOND, BIT,
MetaBAT, external FeGenie clones, and additional model downloads are not needed.
See `environment.yml`, `environment-plots.yml`, and `pyproject.toml` for bounds.

Omit `--plots` for tables only. Existing environments are not deleted; choose
another `--name`, or use `--install-only` after provisioning dependencies.
The installer does not alter global channel settings or shell startup files.
`--dry-run` shows the selected specification. Normal package updates require
reinstalling; use `pip install -e .` only when intentionally developing in place.

Libraries resolve from the checkout or installed data package. An explicit
`-hmm_dir` takes precedence over the library's optional `*_hmms` environment
variable; the bundled library is the default. Remove stale overrides from older
installations. All built-in overrides must match the versioned manifests.

## Three input modes, consistent across every genie

```bash
# Protein FASTA: retain original IDs; no genomic clustering.
python MagicLamp.py GramGenie -bin_dir proteins -bin_ext faa --orfs -out gram_proteins

# Annotated GenBank: use CDS translations and cluster by genomic gene order.
python MagicLamp.py GramGenie -bin_dir annotations -bin_ext gbk --gbk -out gram_genbank

# Nucleotide FASTA: predict proteins with Prodigal, then cluster.
python MagicLamp.py FeGenie -bin_dir genomes -bin_ext fna -out iron_genomes -t 4
```

`--orfs` and `--gbk` are mutually exclusive. `--proteins` is an alias for
`--orfs`. Protein names never imply genomic positions, including Prodigal-like
names. GenBank input uses existing `/translation` qualifiers; pseudogenes and
untranslated CDSs are skipped but still contribute to gene-order gaps. It does
not call Prodigal or translate missing CDSs automatically.

Clustering uses sample and contig boundaries and the genomic CDS order, with
`-d 5` allowing up to five intervening CDSs between neighboring hits. Clusters
are connected neighborhoods, can include singletons, ignore strand, and do not
wrap across a circular chromosome origin. Use `--no-cluster` to disable them for
GenBank/nucleotide input. Protein input always disables them.

Use `--meta` only for Prodigal metagenomic mode on nucleotide FASTA. Other common
options are `-t`, `-cat`, `--norm`, `-inflation` (default 1000 proteins), `-ref`
(reference FASTA for separate BLASTP evidence), and `--makeplots`. Run
`python MagicLamp.py GENIE --help` for the exact options. Output directories must
be new or empty. External command failures stop the run and record failure;
unknown options are rejected.

## Maintained commands and models

| Command | Library | Profiles | Threshold policy |
|---|---|---:|---|
| FeGenie | `hmms/iron` | 192 | Curated sequence score + E ≤ 0.1 |
| LithoGenie | `hmms/litho` | 138 | Curated sequence score + E ≤ 0.1 |
| Lucifer | `hmms/lux` | 19 | Curated sequence score + E ≤ 0.1 |
| ATPGenie | `hmms/atp` | 17 + 4 classifiers | TC; separate proteolipid classification |
| PortGenie | `hmms/portna` | 21 | TC |
| RnfGenie | `hmms/rnf` | 6 | TC |
| MotiliGenie | `hmms/motility` | 61 | TC |
| RiboGenie | `hmms/ribo` | 266 | TC |
| ResistiGenie | `hmms/resist` | 55 | TC |
| SporeGenie | `hmms/spore` | 83 | TC |
| GasGenie | `hmms/gas` | 26 | Curated sequence score + E ≤ 0.1 |
| PlasticGenie | `hmms/plastics` | 50 | TC |
| GramGenie | `hmms/gram` | 6 | Pfam GA, sequence and domain |
| OmniGenie | Selected supported library | As above | Same engine and policy |
| HmmGenie | User-supplied library | Variable | GA/TC; explicit E-value otherwise |

GA and TC policies apply both sequence and domain thresholds. Score policies
retain their original curated sequence cutoffs, without inventing domain
thresholds. Exact files, scores, categories, annotation sources, and SHA-256
checksums are versioned under `latest/manifests/`. Only manifest-listed models
are active; other source profiles under `hmms/` are not automatically searched.
24 FeGenie files without a matching curated cutoff are listed in
`latest/manifests/FeGenie_excluded.tsv` and are inactive. Missing legacy category
annotations are explicitly `unclassified`, with HMM header descriptions retained.
No new phenotype-specific cutoffs have been calibrated in this release.

```bash
python MagicLamp.py OmniGenie -genie atp -bin_dir proteins -bin_ext faa --orfs -out atp_out
python MagicLamp.py HmmGenie -hmm_dir my_profiles -bin_dir proteins -bin_ext faa --orfs -out custom_out
# If custom profiles have no GA/TC, choose an E-value explicitly:
python MagicLamp.py HmmGenie -hmm_dir my_profiles --cutoff evalue -eval 1e-10 -bin_dir proteins -bin_ext faa --orfs -out custom_evalue_out
```

OmniGenie resolves the same named libraries, including nested ATP/iron layouts;
`-genie antibiotics` and `-genie btex` are rejected. HmmGenie accepts complete,
single-profile, amino-acid HMMER3 `.hmm`/`.HMM` files recursively. Custom libraries
are explicitly user-selected and are not certified marker sets.

FeGenie runs the bundled Geobacter and Thermincola reference FASTAs by default;
`--skip-reference` explicitly skips those searches. `--heme` and `--hematite`
report legacy motif candidates across proteins. ATPGenie retains its two-pass
subunit/proteolipid workflow; RnfGenie reports local six-subunit completeness.
See [migration and scientific interpretation](docs/MIGRATION.md) before comparing
with old results.

## Common outputs

| File | Meaning |
|---|---|
| `summary.csv` | One row per passing sample/protein/model; original ID, annotation, sequence, scores, coordinates/context when available |
| `domains.tsv` | All passing domain alignments, including repeated domains |
| `counts.csv` | Protein counts by model and sample, including zero-hit samples |
| `categories.csv` | Distinct proteins per category and sample |
| `samples.tsv` | Protein totals, hit totals, and sample search status |
| `clusters.tsv` | Coordinate-supported neighborhoods; header only for `--orfs` |
| `reference_hits.tsv` | Separate BLASTP reference evidence; not merged into HMM calls |
| `run.json` | Inputs/models/reference checksums, software versions, commands, parameters, and completion/failure status |
| `ORF_calls/`, `HMM_results/` | Search sequences, ID maps, predictions, raw search output, logs |
| `context_supported.csv` | FeGenie subset with explicit supporting context/sequence labels |
| `motif_hits.tsv` | FeGenie motif candidates, when requested |
| `heatmap.png` | Optional plot of the same `counts.csv` values |

`GENIE-summary.csv`, `GENIE.heatmap.csv`, and `GENIE-run.json` aliases use a
lowercase genie name. Alias filenames are familiar; their schemas follow latest.
Counts represent proteins, not domain repeats. A protein can count in multiple
models/categories. `--norm` reports proteins per `-inflation` total input/predicted
proteins. Zero-hit samples have zeros; samples without translated proteins have
blank normalized values because the denominator is zero.

A passing HMM is evidence of sequence similarity, not proof of drug resistance,
EET, substrate specificity, iron metabolism, or a Gram-positive phenotype.
GramGenie's six requested Pfam profiles are cell-wall-associated domains and
sortase; individual matches do not establish cell-wall anchoring or taxonomy.
ResistiGenie's historical stress/resistance annotations also need biological
interpretation; this release does not convert them into validated AMR calls.

## Development and status

```bash
python MagicLamp.py validate
python -m unittest discover -s tests -v
```

CI installs HMMER, Prodigal, BLASTP, and Biopython and checks all maintained
commands with protein and GenBank inputs. Tests also cover profile parsing,
threshold boundaries, repeat retention, zero samples, coordinate gaps, exclusions,
failed runs, and installer registration. They validate implementation, not
biological sensitivity/specificity. Optional plots are tested when installed.

To refresh GramGenie's six original Pfam models:
`python scripts/update_gram_hmms.py --update`. Review and commit both the library
manifest and `latest/manifests/GramGenie.tsv`; without `--update`, changed upstream
checksums are rejected. Other libraries require a reviewed manifest update.

- [Migration and interpretation](docs/MIGRATION.md)
- [Baseline repository audit](docs/REPOSITORY_AUDIT.md)
- [Consolidation verification](docs/LATEST_VALIDATION.md)
- [Historical material and status](legacy/README.md)

The repository retains its existing [license](LICENSE). Preserve original model
attribution and check upstream terms when redistributing HMM collections.
