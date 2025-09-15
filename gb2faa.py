#!/usr/bin/env python3
import argparse
from Bio import SeqIO


def extract_proteins_or_contigs(genbank_file, output_fasta, output_type):
    has_proteins = False
    records = list(SeqIO.parse(genbank_file, "genbank"))
    out = open(output_type, "w")

    # Determine if any sequence contains proteins before processing
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                has_proteins = True
                break
        if has_proteins:
            break

    with open(output_fasta, 'w') as fasta_out:
        if has_proteins:
            out.write("proteins")
            for record in records:
                counter = 0
                for feature in record.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        counter += 1
                        locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                        # product = feature.qualifiers.get("product", ["unknown"])[0]
                        protein_seq = feature.qualifiers["translation"][0]
                        # fasta_out.write(f">{record.name};{locus_tag};{product};{str(counter)}\n{protein_seq}\n")
                        fasta_out.write(f">{record.name};{locus_tag};{counter}\n{protein_seq}\n")
        else:
            out.write("contigs")
            for record in records:
                fasta_out.write(f">{record.name}\n{record.seq}\n")
        out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein sequences from a GenBank file and save them in FASTA format. If no proteins are found, output contigs instead.")
    parser.add_argument("input", help="Path to the input GenBank file")
    parser.add_argument("output", help="Path to the output FASTA file")
    parser.add_argument("outputType", help="Path to the output type file")

    args = parser.parse_args()
    extract_proteins_or_contigs(args.input, args.output, args.outputType)
