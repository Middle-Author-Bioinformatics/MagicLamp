# MagicLamp repository audit

**Date:** 18 November 2026. **Baseline:**
[`bbe400139f24724cec09ce9314cadecaa18afbc4`](https://github.com/Middle-Author-Bioinformatics/MagicLamp/tree/bbe400139f24724cec09ce9314cadecaa18afbc4).

## Main conclusion

The repository has one active dispatcher, but its organization and successful
process exit codes give a stronger impression of reliability than the evidence
supports. `MagicLamp.py` imports `genies_v3`; the other genie directories are not
used by this entry point. Nevertheless, some active modules have reporting or
portability defects, and some inactive libraries are unusable as distributed.
The priority is to make results complete and trustworthy before reorganizing
directories or presenting everything as a supported release.

The most important findings are:

1. **FeGenie is not portable:** the active implementation overrides configured
   library locations with `/home/ark/...` paths and failed the real-tool smoke test.
2. **Annotation tables and HMM libraries disagree:** most notably AbxGenie has
   59 profiles but only 19 matching metadata entries. A real gyrB hit was written
   as `unknown` and then disappeared from its heatmap.
3. **Zero-hit samples disappear from legacy heatmaps.** An absent sample column
   is not interchangeable with a confirmed zero; this can corrupt cohort comparisons.
4. **OmniGenie can report success without searching anything:** its documented
   ATP example points at a nested library, while the runner reads only flat files.
5. **30 BTEX files fail HMMER parsing.** The included splitter drops the HMMER
   format header; the valid combined file contains multiple profiles.
6. **Setup is destructive and environment-specific:** it removes an existing
   `magiclamp` environment, changes global channel settings, pins Python 3.7, and
   derives activation paths from the caller's working directory.

These findings do not invalidate every existing result. They mean that module,
input mode, library version, and output completeness need to be checked explicitly.

## Scope and evidence

This audit examined tracked source, dispatch, installation, HMM layouts,
embedded thresholds, annotation mappings, duplicate scripts/archives, and
observable command-line behavior at the pinned snapshot. Baseline inventory:

| Item | Count |
|---|---:|
| Tracked files | 2,626 |
| Tracked bytes | 356,194,235 (about 340 MiB) |
| Python files | 62 |
| Lowercase `.hmm` files under `hmms/` | 1,282 |
| Additional lowercase `.hmm` files under `depracated/` | 216 |
| Additional uppercase `.HMM` files | 49 |
| HMM library directories containing lowercase `.hmm` files | 22 |
| Existing genies dispatched before this change | 15 |
| Tracked bytecode files | 46 |

All 62 Python files parsed under Python 3.12; this is a syntax check, not a
runtime guarantee. All 15 pre-existing dispatcher commands produced help when
the expected HMM environment variables and repository PATH were supplied.

Real-tool smoke tests used HMMER **3.4**, Prodigal **2.6.3**, and the first
30,000 nucleotides of *Bacillus subtilis* 168 chromosome **NC_000964.3**, with
`--meta -t 1`. Prodigal predicted 25 proteins. This deliberately small fragment
is a pipeline check, not a test of whole-genome trait completeness. Each run
used its own working/output directory. The existing genie implementations were
not repaired before testing. HmmGenie used the six new Pfam profiles as its
custom library; OmniGenie used the README's original `-genie atp` example.

Inspect the underlying evidence:

- [Every baseline Python file and its dispatch status](audit/genie_inventory.tsv)
- [Every lowercase HMM file, checksum, accession, and threshold fields](audit/hmm_inventory.tsv)
- [Runtime outcomes, input checksums, metadata mismatches, archive comparisons, and limitations](audit/validation.json)

The inventories describe the baseline; GramGenie's new six-model inventory is
separate in [`hmms/gram/models.tsv`](../hmms/gram/models.tsv).

Not tested exhaustively: full conda solves, R figures, DIAMOND/BLAST confirmation
modes, every GenBank mode in the older runners, all platforms, every old standalone
script, or biological sensitivity/specificity. This is a source and runtime audit
of a specific snapshot, not a complete review of the repository's history. No
automatic changes were made to old analyses, old libraries, or the user's environments.

## Which versions are in use?

| Location | Baseline Python files | Actual role | Recommended treatment |
|---|---:|---|---|
| `genies_v3/` | 16, including `__init__.py` | Sole package imported by `MagicLamp.py`; 15 dispatched modules | Canonical development location for now; support status must be assigned per module |
| `genies_v2/` | 13 | Not imported by the dispatcher | Inactive snapshot; preserve until differences and external callers are checked |
| `genies/` | 13 | Not imported by the dispatcher | Inactive older implementation; some modules have no v3 equivalent |
| `depracated/` | 16 | Not imported by the dispatcher | Historical/experimental collection, not automatically safe to delete |

Seven v2 files are byte-identical to their v3 counterparts: **ATPGenie,
FeGenie, HmmGenie, LithoGenie, Lucifer, OmniGenie, and RnfGenie**. Other overlapping
copies differ, so a blanket deletion or "latest folder is better" rule is unsafe.
The baseline's latest commit updates `depracated/CircGenie.py`; that is a concrete
reason not to assume every file there is abandoned. Maintainer intent should be
recorded before declaring those unique modules deprecated.

Legacy-only names include CircGenie, MagnetoGenie, MnGenie, PolGenie, RosGenie,
WspGenie, YfGenie, and the separate HmmGenie alignment/build scripts. Their absence
from the dispatcher is verified; their biological validity and standalone runtime
status were not established by this audit.

## Active module status

"Smoke completed" means the tested command reached output generation. It does
not mean all flags work or that its annotations are scientifically validated.

| Module | HMM location | Observed result / practical status |
|---|---|---|
| FeGenie | `hmms/iron` (nested) | **Blocked on another machine:** hard-coded path caused `FileNotFoundError` |
| LithoGenie | `hmms/litho` | Smoke completed; no-hit sample missing from matrix |
| Lucifer | `hmms/lux` | Smoke completed; no-hit sample missing from matrix |
| ATPGenie | `hmms/atp/{f-type,v-a-type,ion-preference}` | Dedicated runner completed; no-hit sample missing; specialized two-pass scoring |
| PortGenie | `hmms/portna` | Smoke completed; no-hit sample missing; legacy GenBank path bug visible in source |
| RnfGenie | `hmms/rnf` | Smoke completed; no-hit sample missing |
| AbxGenie | `hmms/antibiotics` | **Reporting defect:** detected gyrB profile became `unknown` and was omitted from matrix |
| MotiliGenie | `hmms/motility` | Smoke completed; no-hit sample missing; one library/metadata mismatch |
| RiboGenie | `hmms/ribo` | Smoke completed with a hit; 49 profiles lack metadata coverage |
| ResistiGenie | `hmms/resist` | Smoke completed with two hits; metadata/library mismatches remain |
| SporeGenie | `hmms/spore` | Smoke completed; no-hit sample missing; six profiles lack metadata |
| GasGenie | `hmms/gas` | Smoke completed; no-hit sample missing; flat scan does not recurse into additional sublibraries |
| PlasticGenie | `hmms/plastics` | Smoke completed with empty output; flat scan sees 50 of 103 recursive profile files |
| OmniGenie | Selected library | **Unsupported nested-library example:** ATP run returned 0 with no HMM search results |
| HmmGenie | User-provided | Smoke completed with a LysM hit; final summary has no header and contains `EMPTY` locus mapping |
| **GramGenie (new)** | `hmms/gram` | All six synthetic positive controls and real-tool integration tests pass; biological benchmarking still needed |

The baseline has no tests/CI configuration or published GitHub releases. No
Git tags were returned by `git ls-remote --tags`. Four branches were visible
(`master`, `patch-1`, `patch-2`, `patch-3`) and no open pull requests at audit start.
This change adds focused CI for GramGenie, not suite-wide certification.

## Findings and fixes

### P1 — FeGenie ignores its configured HMM location

In [`genies_v3/FeGenie.py`](../genies_v3/FeGenie.py), the environment-derived
directory is read, but the cutoff file at baseline line 376 and library root
at line 569 are reset to `/home/ark/MAB/bin/MagicLamp/hmms/iron...`. The test failed
opening the first of those paths, even with a valid `iron_hmms` value.

**Fix:** resolve the root from an explicit argument, then environment, then the
installed package. Derive the cutoff file from that root. Test from a working
directory outside the checkout and on a machine without `/home/ark`.

### P1 — Library/annotation mapping drift causes lost or misleading results

Comparing profile accessions to each active runner's `HMM_META` gives:

| Runner | Profiles in library | Profiles without metadata | Metadata IDs without a profile |
|---|---:|---:|---:|
| AbxGenie | 59 | **40** | 0 |
| RiboGenie | 266 | **49** | 0 |
| ResistiGenie | 55 | **6** | **33** |
| SporeGenie | 83 | **6** | 0 |
| MotiliGenie | 61 | **1** | **1** |
| PortGenie | 21 | 0 | 0 |
| RnfGenie | 6 | 0 | 0 |

The AbxGenie test provides a reproducible example: `TIGR01059.1.hmm` has
`DESC JCVI: DNA gyrase subunit B`; it matched at 1075.6 bits. The summary assigned
unknown symbol/description/category because the accession is not in `HMM_META`.
The heatmap builder then skipped it, leaving no sample column despite a detected
hit. Full missing-accession lists are in `audit/validation.json`.

There is also a scientific interpretation problem: **DNA gyrase is an antibiotic
target, not a standalone antibiotic-biosynthesis marker.** ResistiGenie and
SporeGenie likewise need evidence-level categories so that housekeeping genes,
stress-response proteins, mutation-dependent targets, and specific acquired
resistance determinants are not treated equivalently. Presence of a domain or
general pathway component is not proof of the phenotype named by a genie.

**Fix:** move model-to-annotation mappings into versioned manifests. Fail
validation when a loaded model lacks a deliberate annotation. Either remove
out-of-scope models after scientific review or label them as contextual evidence.
Never silently discard a detected profile during aggregation.

### P1 — Heatmaps lose zero-hit samples

Many legacy runners derive sample columns from dictionaries populated only by
hits. In this smoke test, 10 dedicated existing modules produced heatmaps with
no named input sample column, including AbxGenie, where an unclassified hit was
also present. The generic ATP route through OmniGenie did so as well.

**Fix:** initialize the full input sample list before searching; write explicit
zero counts for searched samples with no hits, and separate failed/unsearched
samples from biological zeros. Add tests for all-zero cohorts, mixed positive/
negative cohorts, no ORFs, missing executables, and zero normalization denominators.

### P1 — External command failures and invalid requests can look successful

Older modules use unchecked `os.system`, fixed working-directory files such as
`HMMlib.txt`/`mainDir.txt`, shell-interpolated paths, and implicit reuse of existing
ORF files. These mechanisms risk space-sensitive failures, concurrent-run races,
and reuse of stale ORFs when inputs change. Many validation paths use bare
`SystemExit`, giving a successful exit code. `parse_known_args` silently accepts
unknown options. Most older modules also retain only a best HMM per ORF, which
is inappropriate for domain-architecture/repeat analyses.

**Fix:** use argument-list `subprocess.run(check=True)`, isolate intermediates
under each run directory, return nonzero on failure, reject unsupported options,
and use input/model checksums for explicit resume. Keep multiple nonredundant
domains when the scientific question requires them. GramGenie implements these
behaviors without changing old output semantics.

### P1 — Installation can remove an existing environment

[`setup.sh`](../setup.sh) removes `magiclamp` before solving and again on fallback,
mutates global conda channels, pins Python 3.7, and lacks a standalone environment
specification/lock. Its success message overstates runtime readiness.
[`magiclamp.paths`](../magiclamp.paths) uses `pwd` rather than its own location,
so calling setup from another directory writes incorrect HMM roots. It writes
a generic `activate.d/env_vars.sh` and has no matching deactivation restoration.

**Fix:** provide versioned core and optional plotting environment YAMLs, select
a tested current Python, avoid environment deletion by default, and make replacement
explicit. Resolve paths from the script. Use a MagicLamp-specific activation file
with a paired deactivation hook. Keep plotting optional and test the minimal setup
independently. The audit did not run the destructive installer.

### P2 — Generic runner support and directory layout disagree

OmniGenie lists only immediate `.hmm` files. The ATP root has none; all 21 profiles
are nested. Its original README example therefore returned success without a
search. Iron is also nested. The fallback is gated by an `rscripts.txt` file that
the script itself just created, so an empty HMM environment variable is not a
reliable trigger for fallback. Correctly configured flat libraries are a separate
case and are not disproved by the ATP failure.

**Fix:** either restrict generic mode to validated flat libraries or give it a
manifest describing exact files and each profile's cutoff policy. A universal
`--cut_tc` setting cannot safely replace every genie's specialized scoring. The
README now stops promising nested-library support and uses a flat Rnf example.

### P2 — Broken BTEX splitter and ambiguous multi-profile files

`hmmstat` rejected 30 split files in `hmms/btex`: their first tag is `NAME`, not
`HMMER3/...`. [`split_hmms.sh`](../hmms/btex/split_hmms.sh) starts writing only after
seeing `NAME`, losing the preceding format header. `all_models.hmm` parses as a
multi-profile library. Eight CRISPR `.hmm` files also contain multiple models and
parse successfully; multi-model files are not inherently invalid, but runners
must not assume one file equals one profile. Uppercase `.HMM` files are also
ignored by runners using case-sensitive `.endswith('.hmm')`.

**Fix:** regenerate split profiles from a valid original using `hmmfetch` or a
record-aware splitter; require `hmmstat` validation in CI. Enumerate actual model
records, canonicalize extensions, and distinguish source material from runnable
libraries. Do not repair missing headers by guesswork without validating records.

### P2 — Legacy GenBank support is inconsistent

Several active runners, including PortGenie, invoke `gtt-genbank-to-AA-seqs` with
an output path under the output directory but then open that filename under the
input directory. The helper is not explicitly installed by the setup's core
package list. The separate `gb2faa.py` helper imports Biopython but has an absolute
`/home/ark/.../python` shebang, while Biopython is not an explicit core dependency.

**Fix:** use a shared, tested protein-extraction function with explicit dependency
checks and consistent paths/IDs; test compound CDS features, record boundaries,
pseudogenes, and absent translations. GramGenie's optional GenBank mode uses
Biopython directly and was separately exercised with a translated CDS fixture.

### P2 — Repository packaging and documentation obscure support status

The checkout contains three ZIP bundles, inactive script trees, 46 tracked
bytecode files, 18 `.DS_Store` files, six IDE files, and four FUSE artifacts.
`genies.zip` and `genies_v2.zip` duplicate 21 and 34 on-disk files respectively.
`hmms.zip` contains 2,189 non-directory/non-`__MACOSX` entries: 2,180 match files
at their checkout paths, one differs, and eight are not present at those paths.
It is not a reliable replacement for the current library tree.

The baseline README omitted dispatched GasGenie/PlasticGenie, pointed installation
at a different fork, treated some flags as universally supported, and lacked a
per-version compatibility inventory. There was no suite-wide automated test
gate, published release, or model-refresh policy.

**Fix:** stop tracking generated/cache files and expand `.gitignore`; move
release bundles out of the source tree; document one canonical runtime package.
Compare and archive unique older implementations before removing duplicates.
Use a versioned catalogue with support status, inputs, outputs, dependencies,
model sources, cutoff policy, and last passing test. Tag actual tested releases.

## Recommended cleanup sequence

| Priority | Work | Completion criterion |
|---|---|---|
| 1 — Correctness | Repair FeGenie paths; reconcile model annotations; preserve zero-hit columns; fail on missing tools/libraries | Positive, negative, empty, and failing-run fixtures pass for each supported module; no unclassified model is silently dropped |
| 2 — Installation | Non-destructive environment specifications; script-relative paths; dependency preflight | Fresh installation works from outside the checkout on Linux; optional plots are independently tested |
| 3 — Runtime consistency | Shared FASTA/GenBank, HMMER, provenance, and CSV helpers; gradual migration | Stable per-module outputs with regression tests; domain and full-protein workflows remain scientifically distinct |
| 4 — Organization | Canonical package; reviewed legacy archive; remove cache/duplicate ZIPs | Every tracked runtime/library file has an owner, status, and purpose; unique old functionality is preserved |
| 5 — Release management | Catalogue, supported-version policy, CI matrix, tagged release, changelog | A release is reproducible and its support claims match passing checks |

Keep `genies_v3` as the canonical package during repairs. A later rename to a
stable package name such as `magiclamp/genies` can happen as a dedicated migration
with compatibility imports. Do not combine that move with broad algorithm changes.
The desired structure is one runtime package, one manifest per HMM library,
shared utilities, tested fixtures, documentation, and a clearly labeled legacy
archive; generated binaries and release ZIPs belong outside the source tree.

## What this change implements

- A dedicated **GramGenie** with the six requested original Pfam profiles,
  sequence/domain gathering cutoffs, repeat and multidomain retention, input/
  model checksums, explicit failures, and zero-hit sample preservation.
- Dispatcher help, activation-path registration, detailed usage documentation,
  and a checked model-refresh script.
- Focused unit and real-HMMER/Prodigal integration tests and a CI workflow.
- This audit, file/model inventories, and corrected README statements about
  active versions, generic-library support, and installation readiness.

It does **not** repair or delete the older genies, silently reclassify their
scientific markers, or claim that a smoke test establishes trait prediction
accuracy. Those are explicitly prioritized follow-up changes above.
