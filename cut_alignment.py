import argparse
from Bio import AlignIO, SeqIO

def main(alignment_in, alignment_out, start, stop, format):
    if start:
        start = start-1     #Python indexes start at 0 instead of 1
    alignment = AlignIO.read(alignment_in, format)
    cut_alignment = alignment[:, start:stop]
    SeqIO.write(cut_alignment, alignment_out, format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_in", type=str, help="Input alignment")
    parser.add_argument("alignment_out", type=str, help="Output alignment")
    parser.add_argument("--start", type=int, help="First column to include in cut. Default = 1 (first index of the alignment)", default=None)
    parser.add_argument("--stop", type=int, help="Last column to include in cut. Default = End of the alignment", default=None)
    parser.add_argument("--format", type=str, default="fasta", help="Format of the alignment (Default: fasta)")
    args = parser.parse_args()
    main(args.alignment_in, args.alignment_out, args.start, args.stop, args.format)