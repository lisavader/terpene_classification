import argparse
import glob
import ast

from InterproRecord import InterproRecord
from fasta_parsing import write_fasta

def select_records(records, selection_dict):
    selected_records=[]
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
            selected_records.append(record)
    return selected_records

def main(json_dir, selection, fasta_out):
    json_files = glob.glob(json_dir+"/*.json")
    records=[]
    for json_file in json_files:
        record = InterproRecord.from_json(json_file)
        records.append(record)
    with open(selection, "r") as selection_data:
        selection_dict = ast.literal_eval(selection_data.read())
    selected_records = select_records(records, selection_dict)
    selected_accessions = []
    selected_seqs = []
    for selected_record in selected_records:
        selected_accessions.append(selected_record.accession)
        selected_seqs.append(selected_record.seq)
    write_fasta(selected_accessions, selected_seqs, fasta_out)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("json_dir", type=str, help="Path to directory with json files")
    parser.add_argument("selection", type=str, help="Path to file with selection criteria in dict format")
    parser.add_argument("fasta_out", type=str, help="Path to write fasta output file")
    args = parser.parse_args()
    #Run the main script
    main(args.json_dir, args.selection, args.fasta_out)

