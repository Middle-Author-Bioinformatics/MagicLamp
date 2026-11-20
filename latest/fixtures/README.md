# Offline installation fixture

`B_subtilis_first_30kb.fna` contains positions 1–30000 (1-based, inclusive)
of the Bacillus subtilis 168 reference chromosome NC_000964.3. It is a small
pipeline test, not a complete genome or a biological benchmark. It is the same
fragment used in `docs/audit/latest_nucleotide_smoke.json`.

`magiclamp doctor --smoke-test` runs Prodigal and the six GramGenie profiles
on this fragment, checks the expected 25 predicted proteins and at least one
passing GramGenie hit, then checks protein-only and GenBank input separately.
The fixture is bundled in the installation; no download is needed at runtime.
