import argparse
import Bio.Align
from Bio import Phylo, AlignIO

def find_duplicates(tree_in, tree_format, cut_off):
    duplicates = []
    tree=Phylo.read(tree_in, tree_format)
    leaves = tree.get_terminals()
    leaf_pairs = list(zip(leaves, leaves[1:]))
    for pair in leaf_pairs:
        if tree.distance(pair[0],pair[1]) <= float(cut_off):
            duplicates.append(pair[0].name)
    return duplicates

def deduplicate_alignment(alignment_in, aln_format, duplicates):
    deduplicated_records=[]
    alignment = AlignIO.read(alignment_in, aln_format)
    for record in alignment:
        if record.id.startswith(tuple(duplicates)):
            pass
        else:
            deduplicated_records.append(record)
    new_alignment = Bio.Align.MultipleSeqAlignment((deduplicated_records))
    return new_alignment

def main(alignment_in, tree_in, alignment_out, cut_off, aln_format, tree_format):
    duplicates = find_duplicates(tree_in, tree_format, cut_off)
    new_alignment = deduplicate_alignment(alignment_in, aln_format, duplicates)
    AlignIO.write(new_alignment, alignment_out, aln_format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_in", type=str, help="Input alignment")
    parser.add_argument("alignment_out", type=str, help="Deduplicated output alignment")
    parser.add_argument("tree_in", type=str, help="Phylogenetic tree from alignment, in newick format")
    parser.add_argument("cut_off", type=str, help="Branch length distance cut-off; if the distance between adjacent leaf nodes is greater or equal to this value,"+
                        "one of the two sequences will be removed.")
    parser.add_argument("--aln_format", type=str, default="fasta", help="Format of the alignment (Default: fasta)")
    parser.add_argument("--tree_format", type=str, default="newick", help="Format of the phylogenetic tree (Default: newick)")
    args = parser.parse_args()
    main(args.alignment_in, args.tree_in, args.alignment_out, args.cut_off, args.aln_format, args.tree_format)