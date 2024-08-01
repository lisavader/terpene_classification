import argparse
import json
import re
import os

from fasta_parsing import read_fasta, write_fasta

def main(fasta_in, fasta_out, json_dir, select_any, select_all, avoid, id_only, deduplicate):
    selected_accessions = []
    fasta_dict = read_fasta(fasta_in)
    for header in fasta_dict.keys():
        header = re.sub("tr\||sp\|", "", header)
        if "WP_" in header or "XP_" in header:
            accession = '_'.join(re.split(r'[_,()| ]',header)[:2])
        else:
            accession = re.split(r'[_,()| ]',header)[0]
        json_file = os.path.join(json_dir, accession+".json")

        if deduplicate:
            select = True
        else:
            select = None
            for i in range(2):
                try:
                    with open(json_file) as f:
                        json_object = json.load(f)
                except FileNotFoundError:
                    base_accession = accession[2:]
                    json_file = os.path.join(json_dir, base_accession+".json")
                    continue
                break

            if select_any:
                values = extract_values(select_any)
                if any(value in str(json_object) for value in values):
                    select = True
                else:
                    select = False
            if select_all:
                values = extract_values(select_all)
                if all(value in str(json_object) for value in values):
                    select = True
                else:
                    select = False
            if avoid:
                values = extract_values(avoid)
                if any(value in str(json_object) for value in values):
                    select = False
                elif select == None:
                    select = True
        if select:
            selected_accessions.append(accession)
    selected_headers = []
    selected_seqs = []
    for header,seq in fasta_dict.items():
        if any(accession in header for accession in selected_accessions):
            if seq not in selected_seqs:
                selected_headers.append(header)
                selected_seqs.append(seq)
            else:
                print("Skipped duplicate sequence: ",header)
    print("Selected "+str(len(selected_headers))+" sequences.\nWriting to "+fasta_out+".")
    write_fasta(selected_headers, selected_seqs, fasta_out, id_only)

def extract_values(argument):
    try:
        with open(argument,"r") as file:
            values = []
            for line in file:
                values.append(line.rstrip())
    except FileNotFoundError:
        values = argument.split(',')
    return values

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to fasta input file")
    parser.add_argument("fasta_out", type=str, help="Path to fasta output file")
    parser.add_argument("json_dir", type=str, help="Path to directory with metadata files in .json format")
    parser.add_argument("--select_any", type=str, help="Values to select for (any has to be present). Provide either a file with one value per line, or a comma separated string of values.")
    parser.add_argument("--select_all", type=str, help="Values to select for (all have to be present). Provide either a file with one value per line, or a comma separated string of values.")
    parser.add_argument("--avoid", type=str, help="Values to avoid (all have to be avoided). Provide either a file with one value per line, or a comma separated string of values." )
    parser.add_argument("--id_only", action='store_true', help="If true: Print only the id in the header")
    parser.add_argument("--deduplicate", action='store_true', help="Don't perform any specific selection, just remove duplicates.")
    args = parser.parse_args()
    #Run the main script
    main(**vars(args))

