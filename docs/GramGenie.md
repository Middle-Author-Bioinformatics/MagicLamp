# GramGenie

GramGenie searches the six requested Pfam models for cell-wall binding,
cell-wall-associated domains, and sortase. It is available through the active
`genies_v3` dispatcher and does not require R, DIAMOND, or a network connection.

| Pfam | Model | Interpretation |
|---|---|---|
| PF01473 | Choline_bind_1 | Putative cell-wall-binding repeat |
| PF04122 | CW_binding_2 | Cell-wall-binding domain 2 |
| PF00395 | SLH | S-layer homology domain |
| PF04203 | Sortase | Anchoring enzyme, not the substrate's anchoring signal |
| PF01476 | LysM | Polysaccharide-binding domain associated with peptidoglycan/chitin |
| PF07501 | G5 | Domain found in cell-surface and other proteins; contextual evidence |

These are domain annotations, not a Gram-stain classifier. A hit alone does not
establish surface localization, covalent cell-wall attachment, or EET activity.
For example, sortase detection does not identify all of its substrates. Signal
peptides, transmembrane helices, LPXTG sorting signals, genomic neighborhoods,
and MHC co-occurrence are separate analyses and are not inferred here.

## Run

From the repository, with HMMER and (for nucleotide inputs) Prodigal on PATH:

```bash
# Genomes / contigs: call proteins with Prodigal, then search.
python MagicLamp.py GramGenie -bin_dir genomes -bin_ext fna -out gram_out -t 8

# Metagenomic contigs / short assemblies.
python MagicLamp.py GramGenie -bin_dir mags -bin_ext fa --meta -out gram_mags -t 8

# Existing proteins: use their IDs directly; do not call Prodigal.
python MagicLamp.py GramGenie -bin_dir proteins -bin_ext faa --orfs -out gram_proteins -t 8

# Normalize heatmap counts to 1,000 input/predicted proteins per sample.
python MagicLamp.py GramGenie -bin_dir proteins -bin_ext faa --orfs --norm -out gram_norm

# Annotated GenBank CDS translations (optional Biopython dependency).
python MagicLamp.py GramGenie -bin_dir annotations -bin_ext gbk --gbk -out gram_gbk
```

`--proteins` is an alias for `--orfs`. For GenBank input, install Biopython
separately if needed (`conda install -c conda-forge biopython`). Only existing
CDS translations are used; pseudogenes and untranslated CDS features are skipped.
GenBank identifiers are prefixed with the record ID. No gene-neighborhood
clustering is performed. Unknown arguments are rejected rather than ignored.

The output directory must be new or empty. Existing results are never reused
implicitly. Paths containing spaces and FASTA IDs containing pipes are supported;
duplicate IDs within a sample are rejected. Input files are processed in sorted
order. `-t` controls threads within each HMMER search; samples run sequentially.

Library discovery is `-hmm_dir` > `$gram_hmms` > bundled `hmms/gram`, relative to
the installed script rather than the working directory. Each library must contain
exactly the expected profile in each of the six `PFxxxxx.hmm` files, including
both GA cutoffs. `magiclamp.paths` exports `gram_hmms` during normal setup. An
existing installation can use the script-relative fallback immediately after
pulling these changes, or refresh its activation hook from the repository root.

## Thresholds and counting

The original Pfam models are bundled unchanged. `hmmsearch --cut_ga` applies
the per-model sequence and domain gathering cutoffs. The parser also checks
both scores explicitly. This matters for short repeats: a protein can require
multiple convincing domains to pass its sequence threshold. No invented common
E-value or bit-score threshold replaces Pfam's curated cutoffs.

All passing domains are retained. Repeated copies of the same domain are
reported separately in the domain table but count once per protein/family in
the summary and heatmap. Different families on the same protein are retained;
there is no global "best HMM per protein" rule. Counts across families therefore
must not be added to obtain the number of unique hit-containing proteins.
No clan-overlap resolution is applied; overlapping family assignments remain
inspectable using their coordinates.

## Outputs

| File | Contents |
|---|---|
| `gramgenie-summary.csv` | One row per sample/protein/Pfam, full protein sequence, sequence score/E-value, and number of passing domains |
| `gramgenie-domains.tsv` | One row per passing domain, including HMM, alignment, and envelope coordinates; all coordinates are 1-based inclusive |
| `gramgenie.heatmap.csv` | Six Pfam rows by sample; unique hit-containing protein counts, or counts per `-inflation` proteins with `--norm` |
| `gramgenie-samples.tsv` | Protein denominator, unique proteins with any hit, and total passing domains per sample |
| `gramgenie-run.json` | Completion/failure state, input and model checksums, model versions/cutoffs, and executed commands |
| `ORF_calls/` | Proteins searched and, for nucleotide input, Prodigal calls/logs |
| `HMM_results/` | Combined library, raw HMMER domain tables, search output including HMMER version, and logs |

Samples without hits are retained as zero columns. A nucleotide sample with no
predicted proteins has zero raw counts, but blank normalized counts because the
denominator is zero. Empty protein FASTA inputs are rejected. Failed external
commands return a nonzero status and leave a failed run manifest; partial files
must not be interpreted as completed output.

## Models and refresh

[`hmms/gram/models.tsv`](../hmms/gram/models.tsv) records the exact profile
versions, both gathering cutoffs, download URLs, dates, and SHA-256 checksums.
Source family pages are `https://www.ebi.ac.uk/interpro/entry/pfam/<accession>/`.
The models were obtained directly from the documented InterPro Pfam API:

`https://www.ebi.ac.uk/interpro/api/entry/pfam/<accession>/?annotation=hmm`

The bundled models are sufficient to run; no download occurs during analysis.
To verify/retrieve the same current upstream models:

```bash
python scripts/update_gram_hmms.py
# Deliberately accept a changed upstream version; review the resulting diff.
python scripts/update_gram_hmms.py --update
```

The API serves the current upstream models, not an immutable version URL. If
upstream changes, the non-update command refuses the changed checksum. The
committed HMMs are the reproducible snapshot. See the Pfam/EMBL-EBI source terms
at https://www.ebi.ac.uk/about/terms-of-use and API documentation at
https://interpro-documentation.readthedocs.io/en/latest/api.html.

## Verification

```bash
python -m unittest discover -s tests -v
```

Real-tool tests require `hmmsearch`, `hmmemit` (HMMER), and Prodigal on PATH.
They exercise all six profiles using HMM-emitted synthetic positive controls,
repeat counting, multidomain proteins, a no-hit sample, normalization, path
handling, and zero-predicted-protein nucleotide input. Such controls validate
the implementation, not biological sensitivity/specificity or localization.
Tests without available optional executables are explicitly skipped.
