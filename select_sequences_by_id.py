import argparse
from fasta_parsing import read_fasta, write_fasta

def main(fasta_in, fasta_out, accessions, id_only):
    with open(accessions,"r") as file:
        accessions = []
        for line in file:
            accessions.append(line.rstrip())
    selected_headers = []
    selected_seqs = []
    fasta_dict = read_fasta(fasta_in)
    for header,seq in fasta_dict.items():
        if header.startswith(tuple(accessions)):
            if seq not in selected_seqs:
                selected_headers.append(header)
                selected_seqs.append(seq)
    write_fasta(selected_headers, selected_seqs, fasta_out, id_only)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to fasta input file")
    parser.add_argument("fasta_out", type=str, help="Path to fasta output file")
    parser.add_argument("accessions", type=str, help="File with accessions to select, separated by new line")
    parser.add_argument("--id_only", action='store_true', help="If true: Print only the id in the header")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in, args.fasta_out, args.accessions, args.id_only)