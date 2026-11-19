# Migrating to latest (4.0.0)

`MagicLamp.py` now imports only `latest/`. There is one shared implementation of
parsing, search, clustering, outputs, and provenance. Per-genie modules are thin
entry points (`python -m latest.GramGenie`, or the root dispatcher). Historical
packages have moved to `legacy/` and are unsupported. AbxGenie and BTEXGenie code,
libraries, and ZIP snapshots are historical material only; neither is registered
or included by OmniGenie. The remaining built-in libraries use explicit manifests.

## Changes affecting results

This is an output/behavior migration, not a promise of identical legacy results.
The bundled model bytes and their selected cutoffs are preserved. Important
intentional differences are:

- **Every passing family is retained.** The legacy best-hit-only dereplication
  could discard independent domain/family annotations. `summary.csv` has one row
  per protein/model and a domain count; `domains.tsv` retains domain repeats.
  The same protein can therefore contribute to multiple models/categories.
- **Counts represent all passing HMM evidence.** Heatmaps no longer silently
  omit unknown categories or zero-hit samples. They do not silently apply
  downstream FeGenie context exclusions. Inspect context columns and subsets
  separately when evaluating pathway potential.
- **Protein FASTA cannot establish neighborhoods.** `--orfs` does not parse
  numeric suffixes as gene order. It produces no clusters or genomic coordinates;
  neighborhood-dependent evidence is `not_assessed`.
- **GenBank uses real feature positions.** Only translated, non-pseudogene CDSs
  are searched; all CDSs contribute to order/gaps. Different records and samples
  cannot cluster together. Coordinates are 1-based inclusive envelopes; a joined
  CDS uses its envelope. Strand does not restrict clustering; circular origin
  wrapping is not implemented. `-d` is now explicitly the number of intervening
  CDSs, rather than an ambiguous difference between ORF suffixes.
- **Categories come from reviewed manifests.** Existing dictionaries/score
  tables provide available metadata; unmatched annotations are `unclassified`
  with the original HMM name/description. These are data gaps, not validated
  functional classifications. GasGenie resolves complete dotted filenames
  against its score table, fixing the former filename/key mismatch.
- **All normalization uses the same denominator and factor.** `--norm`
  defaults to per 1000 proteins. Empty translated samples have undefined (blank)
  normalized counts. Plot values and count-table values agree.
- **Profiles use explicit thresholds and checksums.** Built-ins reject modified
  files until their manifests are reviewed. FeGenie's 24 unmatched-cutoff files
  remain inactive and are listed individually. Runtime concatenation separates
  HMMER3/b and HMMER3/f records because HMMER cannot mix them in one input stream.
  Temporary query names are unique even when original names/accessions repeat;
  source HMMs are not rewritten.
- **Reference BLASTP evidence is separate.** `-ref` and bundled FeGenie
  reference hits are in `reference_hits.tsv`, with identities/alignment lengths
  and scores. They do not overwrite HMM categories or silently create a validated
  operon assignment. The two FeGenie references can be skipped explicitly.
- **Unknown flags now fail.** Historical analysis/plot/export flags without a
  maintained implementation are not silently accepted. No R runtime is used.
  Old output files and plotting scripts are archived, not compatibility APIs.

## Genie-specific interpretation

FeGenie adds explicit candidate evidence to all raw passing hits. Current rules
cover FLEET (five distinct components), magnetosome (five), each DFE/GACE set
(three), FoxABC (two), FoxEYZ (FoxE required), Cyc1/Cyc2 and CymA/Mtr vicinity,
transport/synthesis component counts, and Mtr/Mto component combinations. Cyc2
length/heme and legacy aromatic-pilus sequence patterns are candidate checks.
The rules are implemented together in `latest/context.py`, rather than repeated
legacy CSV filtering passes. They are **not an exact reproduction of every
legacy filter, category reassignment, or reference-operon matching behavior**.
`context_supported.csv` contains only positively labeled candidates; absence
from this subset does not negate an HMM hit. Isolated regulators and many other
models have no positive rule and remain `neighborhood_only`. Motifs and local
component counts alone do not demonstrate a pathway or EET phenotype.

RnfGenie reports how many of the six distinct subunit accessions occur in each
local hit neighborhood. Protein-only inputs cannot receive local completeness.
ATPGenie searches the 17 TC-qualified subunit models, then the four ion-preference
models only against proteolipid candidates (TIGR01260/PF00137). Classifier models
are not independent positive gene calls or heatmap rows. Comparisons are within
the candidate's complex type; missing alternatives stay unresolved. Score-delta
labels (3/10/20 bits) are inherited heuristics, not calibrated probabilities.

## Installation migration

The installer does not delete environments or edit global channel configuration.
Python 3.11+ replaces the unsupported 3.7 pin. Source-relative library resolution
replaces workstation-specific absolute paths and working-directory assumptions.
Optional `*_hmms` overrides still take precedence and must point to the correct
manifest-compatible library root. An environment-local launcher uses that environment’s Python without changing
PATH or forcing HMM variables. Retire old `env_vars.sh`
hooks in existing environments if they still set obsolete library paths, or use
a fresh environment name.

The baseline audit describes the pre-migration snapshot. Its paths/inventories
are historical; use the runtime manifests and current README for supported state.
