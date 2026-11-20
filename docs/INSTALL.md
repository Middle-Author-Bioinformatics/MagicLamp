# Install MagicLamp and run GramGenie

MagicLamp installs its Python code, versioned HMM manifests, HMM profiles, and
FeGenie reference sequences into the selected environment. The installed
`magiclamp` and `MagicLamp.py` commands do not depend on retaining the cloned
source directory. Bundled analyses run offline after installation. Internet
access is needed to clone the repository and obtain installation dependencies.

## Requirements

Use Linux, macOS, or Linux through WSL2 on Windows. Linux x86_64 is the locally
verified installation platform; the conda environment must be solvable on your
platform. Git and an existing Miniforge/conda installation are prerequisites.

| Dependency | Supported version | Purpose |
|---|---|---|
| Python | 3.11–3.13; installer selects 3.11 | MagicLamp runtime |
| HMMER | ≥3.3, <4 | All HMM searches; `hmmemit` for self-test |
| Prodigal | ≥2.6.3, <3 | Gene prediction for nucleotide FASTA |
| NCBI BLAST+ | ≥2.12, <3 | FeGenie and optional reference protein searches |
| Biopython | ≥1.81, <2 | GenBank input and installation test |
| matplotlib / numpy | ≥3.7, <4 / ≥1.24, <3 | Optional heatmap plots |
| pip / setuptools / wheel | ≥23 / ≥77 / ≥0.41 | Package installation/building |

`environment.yml` includes the complete core suite. `environment-plots.yml`
adds plotting. HMMER, Prodigal, and BLAST+ are compiled external programs: pip
alone does not install them. R, BIT, DIAMOND, MetaBAT, a separate FeGenie checkout,
and additional HMM downloads are not needed. AbxGenie and BTEXGenie are excluded.

## 1. Install Miniforge if conda is not already available

Skip this step if `conda --version` works. Miniforge installation guidance and
platform downloads are maintained at <https://github.com/conda-forge/miniforge>.
On Linux or macOS, download the installer matching the machine:

```bash
curl -fsSL -o Miniforge3.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh
```

Follow the installer prompts and open a new terminal. On a Linux system without
Git, install it using the system package manager (Ubuntu/Debian example:
`sudo apt-get update && sudo apt-get install -y git`). GitHub release installers
are for the local operating system/architecture; use WSL2's Linux terminal on
Windows, not a native Windows environment for the bioinformatics binaries.

## 2. Clone and install

```bash
git clone https://github.com/Middle-Author-Bioinformatics/MagicLamp.git
cd MagicLamp
bash setup.sh --name magiclamp --plots
conda activate magiclamp
```

While installation changes are under review, clone their branch instead:

```bash
git clone --branch codex/standalone-installation \
  https://github.com/Middle-Author-Bioinformatics/MagicLamp.git
```

The installer creates a new environment, installs a normal copy of the Python
package and bundled models, and runs an offline end-to-end check. `--plots` is
optional; omit it for core tables without the heatmap image. Existing
environments are never automatically removed. Choose another `--name` if an
environment called `magiclamp` already exists.

The script can be invoked from any working directory by its absolute path.
It does not change global conda channels, modify shell startup files, or add
HMM/PATH activation hooks. `--dry-run` reports the environment/specification.
`--install-only` installs into an already provisioned environment;
`--register-only` remains an alias for older instructions.

## 3. Verify the installation

```bash
magiclamp --version
magiclamp doctor --smoke-test --plots
```

Omit `--plots` from the check if plotting dependencies were not installed. The
check verifies executables/versions, Python imports, all library checksums and
cutoffs, bundled FeGenie references, and actual GramGenie runs using nucleotide,
protein, and GenBank input. Its temporary files are cleaned up automatically.
A failure returns a nonzero exit code, suitable for a future backend startup
check. Plain `magiclamp doctor` skips the end-to-end controls.

A quick standalone check is to leave the clone and run `magiclamp doctor` from
another directory. Normal installs copy the code and models into the environment;
moving or deleting the clone does not break them. Editable development installs
are the explicit exception.

## 4. Run GramGenie on a folder of genomes

Put one uncompressed genome assembly per FASTA file in one directory:

```text
genomes/
  genome_001.fna
  genome_002.fna
  genome_003.fna
```

Multi-contig assemblies belong in a single file per genome. Files must share the
selected extension. Only matching files directly inside the folder are searched;
subdirectories and gzip files are not searched automatically.

For assembled bacterial genomes:

```bash
magiclamp GramGenie \
  -bin_dir /absolute/path/to/genomes \
  -bin_ext fna \
  -out /absolute/path/to/gramgenie_results \
  -t 24 \
  --makeplots
```

Replace the two paths. Use `-bin_ext fasta` or `-bin_ext fa` if that is your
actual extension. Do not add `--orfs` or `--gbk` for nucleotide FASTA: Prodigal
predicts proteins automatically, and genomic neighborhoods are clustered.
Use `--meta` for metagenomic/short-fragment input. `-t` sets HMMER threads per
sample; genomes are processed sequentially. `--makeplots` needs the plotting
installation; omit it otherwise. The output directory must be new or empty.

If you want normalized counts per 1000 predicted proteins, add `--norm`.
Default counts are raw numbers of proteins matching each model.

For protein FASTA (no clustering):

```bash
magiclamp GramGenie -bin_dir /path/to/proteins -bin_ext faa \
  --orfs -out gram_proteins -t 24
```

For annotated GenBank (clustering from CDS coordinates):

```bash
magiclamp GramGenie -bin_dir /path/to/genbank -bin_ext gbk \
  --gbk -out gram_genbank -t 24
```

GenBank files need CDS `/translation` qualifiers; missing translations are not
predicted automatically. `--no-cluster` disables neighborhoods on coordinate
inputs. `--orfs` always disables clustering.

## 5. Read the output

| Output | Contents |
|---|---|
| `summary.csv` | Passing protein/model matches, scores, annotations, sequences, and available coordinates |
| `domains.tsv` | Domain alignments, including repeats within proteins |
| `counts.csv` | Proteins matching each of the six HMMs per genome, including zeros |
| `categories.csv` | Distinct matching proteins per category and genome |
| `clusters.tsv` | Genomic neighborhoods for nucleotide/GenBank input |
| `samples.tsv` | Per-genome protein totals and search status |
| `heatmap.png` | Optional count heatmap |
| `run.json` | Versions, model/input checksums, commands, parameters, and run status |

The six profiles are PF01473, PF04122, PF00395, PF04203, PF01476, and PF07501.
A protein can match several models and contain repeated domains, so model counts
need not add up to the number of distinct matching proteins. GramGenie screens
cell-wall-associated domains/sortase; it does not assign a Gram phenotype.

## Alternative installation and upgrades

For an environment that already provides the external tools:

```bash
python -m pip install .             # Core Python package and bundled data
# Or include Python plotting dependencies:
python -m pip install '.[plots]'
magiclamp doctor --smoke-test
```

To prepare a named environment manually:

```bash
conda env create -n magiclamp -f environment-plots.yml
conda activate magiclamp
python -m pip install --no-deps --no-build-isolation .
magiclamp doctor --smoke-test --plots
```

For development only, use `python -m pip install -e '.[plots]'`; that install
continues to depend on the checkout. To update a normal install, pull the desired
repository revision and repeat `bash setup.sh --name magiclamp --install-only`
(with `--plots` if needed). If dependency specifications changed, first run:

```bash
conda env update -n magiclamp -f environment-plots.yml
bash setup.sh --name magiclamp --install-only --plots
```

An explicit environment update is never run automatically or pruned. Review the
solver plan. Old `*_hmms` environment variables still override the bundled
libraries; unset stale values if `doctor` reports an external missing/changed
library. Models are not fetched at analysis time. Run records capture resolved
versions; the conda specifications are compatibility bounds, not exact locks.
For a deployment, record the tested environment with
`conda list -n magiclamp --explicit > magiclamp-linux-64.lock.txt` (platform
specific), and pin the repository commit or built wheel used for that deployment.
