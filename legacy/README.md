# Historical material — unsupported

The maintained implementation is `../latest/`, selected by `../MagicLamp.py`.
Nothing in this directory is imported or dispatched by latest.

| Location | Status |
|---|---|
| `genies/`, `genies_v2/`, `genies_v3/` | Archived implementations, retained for comparison and reproducibility |
| `depracated/` | Older source and HMM collections; historical spelling retained |
| `hmms/antibiotics/`, `hmms/btex/` | Excluded AbxGenie / BTEXGenie libraries |
| Other `hmms/` subdirectories | Unsupported/source-only historical libraries |
| `bundles/` | Historical ZIP snapshots; never an installation source |
| `rscripts/`, root helper scripts | Historical utilities and plotting, not used by latest |

Archived scripts may contain absolute paths, obsolete dependencies, incorrect
clustering assumptions, or malformed profiles. Archival retention does not
certify functionality. See `../docs/REPOSITORY_AUDIT.md` for the baseline audit
and `../docs/MIGRATION.md` for the maintained replacement.
