import argparse

from InterproRecord import InterproRecord
from fasta_parsing import read_fasta

def write_json(records, json_out):
    for record in records:
        file_path = json_out+"/"+record.accession+".json"
        json_object = record.to_json()
        with open(file_path,"w") as file:
            file.write(json_object)

def main(input, json_dir):
    records = []
    fasta_dict = read_fasta(input)
    for header,seq in fasta_dict.items():
        record = InterproRecord(header,seq)
        record.assign_enzyme_labels()
        record.assign_taxonomic_labels()
        records.append(record)
    write_json(records, json_dir)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="Path to .fasta file with sequences from InterPro")
    parser.add_argument("json_dir", type=str, help="Path to json output directory")
    args = parser.parse_args()
    #Run the main script
    main(args.input, args.json_dir)