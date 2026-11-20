# Consolidation verification

The maintained dispatcher and package are `MagicLamp.py` and `latest/` (4.0.0).
Legacy sources, unsupported libraries, helpers, R scripts, and ZIP snapshots are
archived under `legacy/`. AbxGenie and BTEXGenie have no maintained registry or
OmniGenie entry. Tracked bytecode, editor state, and generated filesystem files
were removed and ignored.

## Verified locally

- **19 unit/integration tests passed** under Python 3.12, using HMMER 3.4,
  Prodigal 2.6.3, BLAST+ 2.17.0, Biopython, and optional matplotlib/numpy.
- Every one of the **15 maintained commands** ran against protein FASTA and
  annotated GenBank fixtures (30 command/mode combinations). Each used a real
  HMM consensus positive control plus a no-hit sample. Tests confirmed hits,
  explicit zero columns, original IDs, absence of inferred protein-only
  neighborhoods, GenBank coordinates and separate contig neighborhoods.
- HMMER parsed all **944 selected profiles** (940 primary + 4 ATP classifiers),
  in runtime combinations separated by threshold policy and HMMER format.
  Checksums and all embedded GA/TC thresholds matched the committed manifests.
- Threshold boundaries reject below-cutoff sequence/domain scores and retain
  exact-cutoff evidence. Repeated domains and independently matching families
  remain in the results; counts use proteins rather than repeats.
- GenBank tests included pseudogene gaps, repeated locus tags on different
  records, reverse-strand features, and `--no-cluster`. Rnf completeness counts
  distinct local subunits; Fe context checks count distinct components.
- ATP classification was checked with real proteolipid searches and a controlled
  missing-alternative case. Opposite-complex scores cannot win the comparison;
  a missing alternative remains unresolved.
- FeGenie ran outside the checkout with bundled HMMs and both bundled reference
  FASTAs, without workstation-specific paths. Its HMMER3/b and HMMER3/f profiles
  are searched separately. Explicit BLASTP evidence is retained independently.
- Failed searches produce failed provenance and no successful summary; existing
  outputs are not reused. Strict arguments, custom-library threshold requirements,
  modified-model rejection, and excluded aliases were checked.
- Protein-count normalization, zero-protein denominators, and optional plotting
  passed. Installer dry-run paths were checked from a different working directory;
  a conda stub verified failure propagation and executable registration, including
  paths with spaces and refusal to overwrite unrelated executables. The launcher
  does not add activation/deactivation hooks or alter PATH.
- Python parsing, Bash syntax, source whitespace, direct/module entry points,
  and the model-refresh script's CLI were checked.

All 15 commands also completed the same **30 kb Bacillus subtilis NC_000964.3
fragment** used in the baseline audit, with 25 Prodigal-predicted proteins per
run. See [captured counts and versions](audit/latest_nucleotide_smoke.json).
Zero hits on a short fragment are expected and are preserved. The baseline
FeGenie portability failure and OmniGenie nested-ATP no-search failure no longer
occur in this check.

## Limits

Consensus controls establish that models are reachable and outputs are correctly
retained; they do not measure sensitivity, specificity, or phenotype prediction.
This work did not recalibrate HMMs, adjudicate every legacy annotation, benchmark
large cohorts, perform a full conda solve, or verify macOS/Windows installations.
The installer behavior tests use a stub; dependencies were exercised through
separately installed real tools. Migration changes, including FeGenie's explicit
candidate context rules and separate reference evidence, are described in
[MIGRATION.md](MIGRATION.md). Legacy final tables are not output-compatible.

GitHub Actions repeats the shared-suite tests under Python 3.11 and 3.12 using
Ubuntu-packaged HMMER, Prodigal, and BLASTP. Its result is reported on the PR.

## Standalone installation (4.0.1)

A fresh Miniforge installation successfully solved `environment-plots.yml`,
created a new Linux x86_64 environment, built/installed the package, and passed
`magiclamp doctor --smoke-test --plots`. The package contains all 944 manifest
profiles, both FeGenie reference FASTAs, and the offline genome fixture. There
are no legacy, AbxGenie, or BTEXGenie libraries in the wheel.

Both installed command names and all three GramGenie input-mode controls passed
with the installation source directory renamed out of the way and with no
source path on PYTHONPATH. Installed code/model paths resolved within the conda
environment. A separate installed FeGenie nucleotide/reference run completed;
`pip check` reported no broken requirements. The expanded 24-test suite passed.
Package and executable versions are captured in
[audit/standalone_installation.json](audit/standalone_installation.json).

These checks supersede the earlier installer/full-conda-solve limitation for
this Linux installation path. macOS/WSL2 installation and biological performance
remain untested. CI additionally builds the source distribution and wheel,
installs the wheel into a separate environment, and runs its offline self-test
outside the checkout on Python 3.11 and 3.12.
