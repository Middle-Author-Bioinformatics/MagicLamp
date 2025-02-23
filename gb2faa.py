#!/usr/bin/env python3
import argparse
from Bio import SeqIO


def extract_proteins_or_contigs(genbank_file, output_fasta, output_type):
    has_proteins = False
    out = open(output_type, "w")
    counter = 0
    with open(output_fasta, 'w') as fasta_out:
        for record in SeqIO.parse(genbank_file, "genbank"):
            contig_name = record.name  # Contig name from LOCUS line
            contig_seq = str(record.seq)
            counter = 0

            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    has_proteins = True
                    counter += 1
                    locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    protein_seq = feature.qualifiers["translation"][0]
                    fasta_out.write(f">{contig_name};{locus_tag};{counter}\n{protein_seq}\n")

            # If no protein sequences were found, write contigs instead
            if not has_proteins:
                counter += 1
                fasta_out.write(f">{contig_name}\n{contig_seq}\n")
    if counter > 0:
        out.write("contigs")
    else:
        out.write("proteins")
    out.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein sequences from a GenBank file and save them in FASTA format. If no proteins are found, output contigs instead.")
    parser.add_argument("input", help="Path to the input GenBank file")
    parser.add_argument("output", help="Path to the output FASTA file")
    parser.add_argument("outputType", help="Path to the output FASTA file")

    args = parser.parse_args()
    extract_proteins_or_contigs(args.input, args.output, args.outputType)
