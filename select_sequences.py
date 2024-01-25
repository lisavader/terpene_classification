import argparse
import json
import jmespath
import re

from fasta_parsing import read_fasta, write_fasta

def main(fasta_in, fasta_out, json_dir, query, id_only):
    selected_accessions = []
    fasta_dict = read_fasta(fasta_in)
    for header in fasta_dict.keys():
        match = re.search(".*?\.\d",header)   #for ncbi ids
        if match:
            accession = match.group(0)
        else:
            accession = re.split(r'[_,()| ]',header)[0]  #for other ids
        json_file = json_dir+accession+".json"

        with open(json_file) as f:
            json_object = json.load(f)
        found = jmespath.search(query, data=json_object)
        print(json_file, found)
        if found:
            accession = jmespath.search("accession", data=json_object)
            selected_accessions.append(accession)
    selected_headers = []
    selected_seqs = []

    for header,seq in fasta_dict.items():
        if header.startswith(tuple(selected_accessions)):
            if seq not in selected_seqs:
                selected_headers.append(header)
                selected_seqs.append(seq)
    write_fasta(selected_headers, selected_seqs, fasta_out, id_only)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to fasta input file")
    parser.add_argument("fasta_out", type=str, help="Path to fasta output file")
    parser.add_argument("json_dir", type=str, help="Path to directory with metadata files in .json format")
    parser.add_argument("query", type=str, help="A query for searching metadata, may use ==, =!, &&, || and (). "+
                        "For example: review_status=='reviewed' && enzyme_type.contains(@,'phytoene synthase')")
    parser.add_argument("--id_only", action='store_true', help="If true: Print only the id in the header")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in, args.fasta_out, args.json_dir, args.query, args.id_only)

