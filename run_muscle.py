import argparse
import subprocess

def main(fasta_in, msa_out):
    command = ["muscle", "-in", fasta_in, "-out", msa_out]
    subprocess.run(command)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to fasta input file")
    parser.add_argument("msa_out", type=str, help="Path to msa output file")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in, args.msa_out)