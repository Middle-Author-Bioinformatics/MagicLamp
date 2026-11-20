"""One explicit catalogue for CLI dispatch, library selection, and installation."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Source checkouts keep hmms/ alongside latest/. Wheels carry the same data in
# a dedicated package, so an installed command never needs the source checkout.
if (ROOT / "hmms").is_dir():
    HMM_ROOT = ROOT / "hmms"
else:
    import magiclamp_models
    HMM_ROOT = Path(magiclamp_models.__file__).resolve().parent
GENIES = {
    "FeGenie": ("iron", "iron_hmms", "Iron genes, heme motifs, and genomic context"),
    "LithoGenie": ("litho", "litho_hmms", "Chemolithotrophic gene profiles"),
    "Lucifer": ("lux", "lux_hmms", "Light sensing and production profiles"),
    "ATPGenie": ("atp", "atp_hmms", "ATP synthase subunits and ion-preference candidates"),
    "PortGenie": ("portna", "portna_hmms", "Sodium-coupled transporter profiles"),
    "RnfGenie": ("rnf", "rnf_hmms", "Rnf subunits and local complex completeness"),
    "MotiliGenie": ("motility", "motility_hmms", "Motility, chemotaxis, and pili profiles"),
    "RiboGenie": ("ribo", "ribo_hmms", "Ribosomal and translation profiles"),
    "ResistiGenie": ("resist", "resist_hmms", "Resistance/stress-associated profile evidence"),
    "SporeGenie": ("spore", "spore_hmms", "Sporulation-associated profiles"),
    "GasGenie": ("gas", "gas_hmms", "Gas-vesicle profiles"),
    "PlasticGenie": ("plastics", "plastic_hmms", "Plastic/pollutant-associated profiles"),
    "GramGenie": ("gram", "gram_hmms", "Cell-wall-associated domains and sortase"),
    "OmniGenie": (None, None, "Run a supported library with the shared engine"),
    "HmmGenie": (None, None, "Run an explicitly supplied HMM library"),
}
EXCLUDED = {"abx", "abxgenie", "antibiotics", "btex", "btexgenie"}


def resolve_name(value):
    lower = value.lower()
    if lower in EXCLUDED:
        raise ValueError("AbxGenie and BTEXGenie are excluded from latest; see legacy/ for historical material")
    for name, (directory, variable, _) in GENIES.items():
        aliases = {name.lower()}
        if directory:
            aliases.update([directory, variable.removesuffix("_hmms")])
        if lower in aliases:
            return name
    raise ValueError("Unknown genie {!r}; run MagicLamp.py help".format(value))
