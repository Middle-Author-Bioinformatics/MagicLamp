"""Sequence parsing and genomic positions, independent of any genie."""
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Protein:
    identifier: str
    sequence: str
    contig: str = ""
    start: object = ""
    end: object = ""
    strand: object = ""
    order: object = ""


def fasta(path, allow_empty=False):
    records = []
    identifier, description, sequence = None, "", []
    seen = set()

    def finish():
        if identifier is None:
            return
        seq = "".join(sequence).upper().rstrip("*")
        if not seq or set(seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError("Empty or invalid protein sequence {} in {}".format(identifier, path))
        if identifier in seen:
            raise ValueError("Duplicate protein ID {} in {}".format(identifier, path))
        seen.add(identifier)
        records.append((identifier, seq, description))

    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                description = line[1:]
                parts = description.split()
                if not parts:
                    raise ValueError("Empty FASTA identifier in {}".format(path))
                identifier, sequence = parts[0], []
            else:
                if identifier is None:
                    raise ValueError("Sequence before FASTA header in {}".format(path))
                sequence.append(line)
    finish()
    if not records and not allow_empty:
        raise ValueError("No protein sequences in {}".format(path))
    return records


def protein_fasta(path):
    # FASTA names, including numeric suffixes, carry no positional meaning.
    return [Protein(identifier, sequence) for identifier, sequence, _ in fasta(path)]


def prodigal_fasta(path):
    proteins, ranks = [], {}
    for identifier, sequence, description in fasta(path, allow_empty=True):
        fields = description.split(" # ")
        if len(fields) < 4:
            raise ValueError("Missing Prodigal coordinates in {}".format(path))
        contig = identifier.rsplit("_", 1)[0]
        ranks[contig] = ranks.get(contig, 0) + 1
        proteins.append(Protein(identifier, sequence, contig, int(fields[1]), int(fields[2]),
                                int(fields[3]), ranks[contig]))
    return proteins


def genbank(path):
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ValueError("--gbk requires Biopython; install the supplied environment") from exc
    proteins, records, identifiers = [], set(), set()
    with Path(path).open() as handle:
        for record in SeqIO.parse(handle, "genbank"):
            if record.id in records:
                raise ValueError("Duplicate GenBank record ID {} in {}".format(record.id, path))
            records.add(record.id)
            features = sorted((f for f in record.features if f.type == "CDS"),
                              key=lambda f: (int(f.location.start), int(f.location.end)))
            # Rank every CDS, including skipped/untranslated/pseudogene CDS, so gaps are real.
            for order, feature in enumerate(features, 1):
                q = feature.qualifiers
                if "pseudo" in q or "pseudogene" in q or not q.get("translation"):
                    continue
                tag = q.get("locus_tag", q.get("protein_id", ["CDS{}".format(order)]))[0]
                identifier = record.id + "|" + tag
                if identifier in identifiers:
                    raise ValueError("Duplicate CDS identifier {} in {}".format(identifier, path))
                identifiers.add(identifier)
                sequence = "".join(q["translation"][0].split()).upper().rstrip("*")
                if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
                    raise ValueError("Invalid CDS translation {} in {}".format(identifier, path))
                proteins.append(Protein(identifier, sequence, record.id, int(feature.location.start) + 1,
                                        int(feature.location.end), feature.location.strand or "", order))
    if not records:
        raise ValueError("No GenBank records in {}".format(path))
    return proteins


def write_search_fasta(path, proteins):
    lookup = {}
    with Path(path).open("w") as handle:
        for index, protein in enumerate(proteins, 1):
            key = "p{:08d}".format(index)
            lookup[key] = protein
            handle.write(">{}\n{}\n".format(key, protein.sequence))
    return lookup
