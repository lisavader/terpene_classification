import argparse
import csv

from InterproRecord import InterproRecord
from fasta_parsing import read_fasta

def write_json(records, json_out):
    for record in records:
        file_path = json_out+"/"+record.accession+".json"
        json_object = record.to_json()
        with open(file_path,"w") as file:
            file.write(json_object)

def write_metadata_summary(records, table_out):
    with open(table_out, 'w') as handle:
        writer = csv.writer(handle)
        header = ["accession","protein_name","enzyme_class","enzyme_subclass","organism_name","organism_division","review_status"]
        writer.writerow(header)
        for record in records:
            row = []
            for attribute in header:
                value = getattr(record, attribute)
                row.append(value)
            writer.writerow(row)

def main(input, json_out, table_out):
    records = []
    fasta_dict = read_fasta(input)
    for header,seq in fasta_dict.items():
        record = InterproRecord(header,seq)
        record.assign_enzyme_labels()
        record.assign_taxonomic_labels()
        records.append(record)
    write_json(records, json_out)
    write_metadata_summary(records, table_out)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="Path to .fasta file with sequences from InterPro")
    parser.add_argument("json_out", type=str, help="Path to json output directory")
    parser.add_argument("table_out", type=str, help="Path to summary output table in .csv format")
    args = parser.parse_args()
    #Run the main script
    main(args.input, args.json_out, args.table_out)