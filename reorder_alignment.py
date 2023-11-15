import argparse
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

def main(alignment_in, alignment_out, accessions, format):
    with open(accessions,"r") as file:
        accession_list = []
        for line in file:
            accession_list.append(line.rstrip())
    alignment = AlignIO.read(alignment_in, format)
    #split the alignment into two groups, based on accessions
    group1 = MultipleSeqAlignment([])
    group2 = MultipleSeqAlignment([])
    for row in alignment:
        if row.id in accession_list:
            group1.append(row)
        else:
            group2.append(row)
    #merge the groups to produce one alignment again
    reordered_alignment = [*group1, *group2]
    SeqIO.write(reordered_alignment, alignment_out, format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_in", type=str, help="Input alignment")
    parser.add_argument("alignment_out", type=str, help="Output alignment")
    parser.add_argument("accessions", type=str, help="List of accessions that should be grouped, separated by new lines")
    parser.add_argument("--format", type=str, default="fasta", help="Format of the alignment (Default: fasta)")
    args = parser.parse_args()
    main(args.alignment_in, args.alignment_out, args.accessions, args.format)