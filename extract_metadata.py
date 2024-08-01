import argparse
from Records import InterproRecord, UniProtRecord
from fasta_parsing import read_fasta

def write_json(records, json_out):
    for record in records:
        file_path = json_out+"/"+record.accession+".json"
        json_object = record.to_json()
        with open(file_path,"w") as file:
            file.write(json_object)

def main(fasta_in, json_dir, database):
    records = []
    fasta_dict = read_fasta(fasta_in)
    for header in fasta_dict.keys():
        if database == "interpro":
            record = InterproRecord(header)
        elif database == "uniprot":
            record = UniProtRecord(header)
        else:
            exit("Error: Not a valid database: "+database)
        record.assign_enzyme_type()
        record.assign_taxonomic_labels()
        records.append(record)
    write_json(records, json_dir)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to .fasta file with sequences downloaded from InterPro or UniProt")
    parser.add_argument("json_dir", type=str, help="Path to json output directory")
    parser.add_argument("database", type=str, help="Origin database. Choose 'interpro' or 'uniprot'.")
    args = parser.parse_args()
    #Run the main script
    main(**vars(args))