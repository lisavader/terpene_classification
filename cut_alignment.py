import argparse
from Bio import AlignIO

def main(alignment_in, alignment_out, regions, format):
    alignment = AlignIO.read(alignment_in, format)
    fields = regions.split(',')
    alignment_regions = []
    for field in fields:
        if ':' in field:
            start,stop = map(int, field.split(':'))
            start = start-1 #Input index starts at 1 instead of 0
            alignment_region = alignment[:, start:stop]
        else:
            alignment_region = alignment[:, int(field)]
        alignment_regions.append(alignment_region)
    #paste the regions together
    merged_alignment = alignment_regions[0]
    for align in alignment_regions[1:]:
        for i, record in enumerate(align):
            merged_alignment[i].seq+=record.seq
    AlignIO.write(merged_alignment, alignment_out, format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_in", type=str, help="Input alignment")
    parser.add_argument("alignment_out", type=str, help="Output alignment")
    parser.add_argument("regions", type=str, help="Which region(s) to cut. If multiple, separate by a comma. E.g. 1:8,10:12,21.")
    parser.add_argument("--format", type=str, default="fasta", help="Format of the alignment (Default: fasta)")
    args = parser.parse_args()
    main(args.alignment_in, args.alignment_out, args.regions, args.format)