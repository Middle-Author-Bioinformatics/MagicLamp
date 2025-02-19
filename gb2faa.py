import argparse
from Bio import SeqIO


def extract_proteins(genbank_file, output_fasta):
    with open(output_fasta, 'w') as fasta_out:
        for record in SeqIO.parse(genbank_file, "genbank"):
            contig_name = record.name  # Contig name from LOCUS line
            counter = 0
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    counter += 1
                    locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    protein_seq = feature.qualifiers["translation"][0]
                    fasta_out.write(f">{contig_name};{locus_tag};{counter}\n{protein_seq}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein sequences from a GenBank file and save them in FASTA format.")
    parser.add_argument("input", help="Path to the input GenBank file")
    parser.add_argument("output", help="Path to the output FASTA file")

    args = parser.parse_args()
    extract_proteins(args.input, args.output)