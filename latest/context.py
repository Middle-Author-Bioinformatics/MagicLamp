"""Genomic neighborhoods and explicitly labeled genie-specific evidence."""
from collections import defaultdict
import re


def assign_clusters(hits, proteins, enabled, max_intervening=5):
    if not enabled:
        for hit in hits:
            hit["cluster_id"] = ""
        return []
    by_contig = defaultdict(set)
    for hit in hits:
        protein = proteins[hit["target"]]
        if protein.contig and protein.order != "":
            by_contig[protein.contig].add(hit["target"])
    clusters, count = [], 0
    for contig, targets in sorted(by_contig.items()):
        ordered = sorted(targets, key=lambda t: proteins[t].order)
        groups = []
        for target in ordered:
            if not groups or proteins[target].order - proteins[groups[-1][-1]].order - 1 > max_intervening:
                groups.append([])
            groups[-1].append(target)
        for group in groups:
            count += 1
            cluster = "C{:06d}".format(count)
            members = set(group)
            for hit in hits:
                if hit["target"] in members:
                    hit["cluster_id"] = cluster
            clusters.append({"cluster_id": cluster, "contig": contig, "start": min(proteins[t].start for t in group),
                             "end": max(proteins[t].end for t in group), "protein_count": len(group),
                             "proteins": ";".join(proteins[t].identifier for t in group)})
    return clusters


def heme_motifs(sequence):
    return sum(len(re.findall("C" + "." * gap + "CH", sequence)) for gap in [2, 3, 4, 14, 15])


def annotate(name, hits, proteins, enabled):
    grouped = defaultdict(list)
    for hit in hits:
        hit["context_status"] = "not_assessed" if not enabled else "neighborhood_only"
        hit["context_annotation"] = ""
        hit["complex_completeness"] = ""
        if hit.get("cluster_id"):
            grouped[hit["cluster_id"]].append(hit)
        if name == "FeGenie":
            sequence = proteins[hit["target"]].sequence
            hit["heme_c_motifs"] = heme_motifs(sequence)
            hit["heme_b_motifs"] = len(re.findall(r"G.[HR].C[PLAV]G", sequence))
    if not enabled:
        return
    if name == "RnfGenie":
        for group in grouped.values():
            subunits = {h["model"].accession.split(".")[0] for h in group}
            complete = len(subunits & {"TIGR0194" + str(i) for i in range(3, 9)})
            for hit in group:
                hit["complex_completeness"] = "{}/6".format(complete)
                hit["context_status"] = "complete_local_set" if complete == 6 else "partial_local_set"
    if name == "FeGenie":
        for group in grouped.values():
            _iron_context(group, proteins)


def _iron_context(group, proteins):
    # Retain raw profile hits; these rules add evidence labels rather than deleting them.
    def key(hit):
        return hit["model"].path.stem.split(".")[0]
    genes = {key(hit) for hit in group}
    categories = defaultdict(set)
    for hit in group:
        categories[hit["model"].category].add(key(hit))
    fleet = {"EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"}
    mam = {"MamA", "MamB", "MamE", "MamK", "MamP", "MamM", "MamQ", "MamI", "MamL", "MamO"}
    dfe1 = {"DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"}
    dfe2 = {"DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"}
    gace = {"GACE_1843", "GACE_1844", "GACE_1845", "GACE_1846", "GACE_1847"}
    mtr = {"MtrA", "MtoA", "MtrB_TIGR03509", "MtrC_TIGR03507"}
    for hit in group:
        gene, category = key(hit), hit["model"].category
        required, observed = None, None
        if gene in fleet:
            required, observed = 5, len(genes & fleet)
        elif gene in mam:
            required, observed = 5, len(genes & mam)
        elif gene in dfe1:
            required, observed = 3, len(genes & dfe1)
        elif gene in dfe2:
            required, observed = 3, len(genes & dfe2)
        elif gene in gace:
            required, observed = 3, len(genes & gace)
        elif gene in {"FoxA", "FoxB", "FoxC"}:
            required, observed = 2, len(genes & {"FoxA", "FoxB", "FoxC"})
        elif gene in {"FoxE", "FoxY", "FoxZ"}:
            required, observed = 1, int("FoxE" in genes)
        elif gene == "Cyc1":
            required, observed = 1, int(bool(genes & {"Cyc2_repCluster1", "Cyc2_repCluster2", "Cyc2_repCluster3"}))
        elif gene == "CymA":
            required, observed = 1, int(bool(genes & mtr))
        elif category in {"iron_aquisition-siderophore_transport_potential", "iron_aquisition-heme_transport",
                          "iron_aquisition-siderophore_transport"}:
            required, observed = 2, len(set().union(*(v for k, v in categories.items() if k in {
                "iron_aquisition-siderophore_transport_potential", "iron_aquisition-heme_transport", "iron_aquisition-siderophore_transport"})))
        elif category == "iron_aquisition-siderophore_synthesis":
            required, observed = 3, len(categories[category])
        elif category in {"iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"}:
            required, observed = 2, len(categories["iron_aquisition-iron_transport"] | categories["iron_aquisition-heme_oxygenase"])
        if required is not None:
            hit["context_status"] = "supported_neighborhood" if observed >= required else "insufficient_context"
            hit["context_annotation"] = "{} distinct required components; {} observed".format(required, observed)
        if gene in {"Cyc2_repCluster1", "Cyc2_repCluster2", "Cyc2_repCluster3"}:
            supported = len(proteins[hit["target"]].sequence) >= 365 and hit["heme_c_motifs"] > 0
            hit["context_status"] = "sequence_supported_candidate" if supported else "insufficient_sequence_evidence"
        if gene in mtr:
            if {"MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507"} <= genes:
                hit["context_status"], hit["context_annotation"] = "supported_neighborhood", "MtrABC-like iron-reduction candidate"
            elif "MtoA" in genes and "MtrB_TIGR03509" in genes and "MtrC_TIGR03507" not in genes:
                hit["context_status"], hit["context_annotation"] = "supported_neighborhood", "MtoAB-like iron-oxidation candidate"
            elif ("MtrB_TIGR03509" in genes and genes & {"MtrA", "MtoA", "MtrC_TIGR03507"}) or "MtrC_TIGR03507" in genes:
                hit["context_status"], hit["context_annotation"] = "partial_neighborhood", "Incomplete porin/cytochrome candidate set"
            else:
                hit["context_status"] = "insufficient_context"
        if gene == "t4ap":
            seq = proteins[hit["target"]].sequence
            aromatics = set("FYW H".replace(" ", ""))
            fraction = sum(aa in aromatics for aa in seq) / len(seq)
            pattern = re.search(r"[FYWH].{22}[FYWH].{2}[FYWH].{4}[FYWH].{17}[FYWH][FYWH].{5}[FYWH]", seq)
            long_gap = re.search(r"[^FYWH]{36,}", seq)
            hit["context_status"] = "sequence_supported_candidate" if fraction > .097 and pattern and not long_gap else "insufficient_sequence_evidence"
