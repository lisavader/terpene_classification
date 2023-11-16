import argparse
import re
from Bio import AlignIO

def main(alignment_in, alignment_out, format):
    alignment = AlignIO.read(alignment_in, format)
    for record in alignment:
        accession = re.split(r'[_,()| ]',record.id)[0]
        record.id = accession
        record.description = accession
    AlignIO.write(alignment, alignment_out, format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_in", type=str, help="Input alignment")
    parser.add_argument("alignment_out", type=str, help="Output alignment")
    parser.add_argument("--format", type=str, default="fasta", help="Format of the alignment (Default: fasta)")
    args = parser.parse_args()
    main(args.alignment_in, args.alignment_out, args.format)