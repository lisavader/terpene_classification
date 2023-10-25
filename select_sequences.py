import argparse
import glob
import ast

from InterproRecord import InterproRecord
from fasta_parsing import read_fasta, write_fasta

def select_accessions(records, selection_dict):
    selected_accessions=[]
    for record in records:
        for attribute_name,selection_values in selection_dict.items():
            match = False   #whether a match is found for this attribute
            attribute_value = getattr(record, attribute_name)
            if isinstance(attribute_value, list):    #check for all list items
                for item in attribute_value:
                    if item in selection_values:
                        match = True
                        break
            else:
                if attribute_value in selection_values:
                    match = True
            if match == False:  #all attributes in the selection dict must yield a match
                break
        if match == True:
            selected_accessions.append(record.accession)
    return selected_accessions

def main(fasta_in, fasta_out, json_dir, selection):
    json_files = glob.glob(json_dir+"/*.json")
    records=[]
    for json_file in json_files:
        record = InterproRecord.from_json(json_file)
        records.append(record)
    with open(selection, "r") as selection_data:
        selection_dict = ast.literal_eval(selection_data.read())
    selected_accessions = select_accessions(records, selection_dict)

    selected_headers = []
    selected_seqs = []
    fasta_dict = read_fasta(fasta_in)
    for header,seq in fasta_dict.items():
        if header.startswith(tuple(selected_accessions)):
            selected_headers.append(header)
            selected_seqs.append(seq)
    write_fasta(selected_headers, selected_seqs, fasta_out)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to fasta input file")
    parser.add_argument("fasta_out", type=str, help="Path to fasta output file")
    parser.add_argument("json_dir", type=str, help="Path to directory with metadata files in .json format")
    parser.add_argument("selection", type=str, help="Path to file with selection criteria in dict format")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in, args.fasta_out, args.json_dir, args.selection)

