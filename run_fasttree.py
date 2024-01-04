import argparse
import subprocess

def main(msa_in, tree_out):
    command = ["FastTree", "-out", tree_out, msa_in]
    subprocess.run(command)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("msa_in", type=str, help="Path to msa input file")
    parser.add_argument("tree_out", type=str, help="Path to tree output file")
    args = parser.parse_args()
    #Run the main script
    main(args.msa_in, args.tree_out)