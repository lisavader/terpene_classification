import argparse
from fasta_parsing import read_fasta, write_fasta
from Bio import Phylo

def find_duplicates(tree_in, tree_format, cut_off):
    duplicates = []
    tree=Phylo.read(tree_in, tree_format)
    leaves = tree.get_terminals()
    leaf_pairs = list(zip(leaves, leaves[1:]))
    for pair in leaf_pairs:
        if tree.distance(pair[0],pair[1]) <= float(cut_off):
            duplicates.append(pair[0].name)
    return duplicates

def main(fasta_in, fasta_out, tree_in, cut_off, tree_format):
    duplicates = find_duplicates(tree_in, tree_format, cut_off)
    selected_headers = []
    selected_seqs = []
    fasta_dict = read_fasta(fasta_in)
    for header,seq in fasta_dict.items():
        if not header.startswith(tuple(duplicates)):
            if seq not in selected_seqs:
                selected_headers.append(header)
                selected_seqs.append(seq)
    print("Removed "+str(len(duplicates))+" sequences.")
    write_fasta(selected_headers, selected_seqs, fasta_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Input fasta file")
    parser.add_argument("fasta_out", type=str, help="Deduplicated fasta file")
    parser.add_argument("tree_in", type=str, help="Phylogenetic tree from sequences, in newick format")
    parser.add_argument("cut_off", type=str, help="Branch length distance cut-off; if the distance between adjacent leaf nodes is greater or equal to this value,"+
                        "one of the two sequences will be removed.")
    parser.add_argument("--tree_format", type=str, default="newick", help="Format of the phylogenetic tree (Default: newick)")
    args = parser.parse_args()
    main(args.fasta_in, args.fasta_out, args.tree_in, args.cut_off, args.tree_format)