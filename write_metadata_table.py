import argparse
import csv
import glob
from Records import ProteinRecord

def write_metadata_summary(records, table_out):
    with open(table_out, 'w') as handle:
        writer = csv.writer(handle)
        header = ["accession","database","protein_name","enzyme_type","organism_name","organism_division","reviewed"]
        writer.writerow(header)
        for record in records:
            row = []
            for attribute in header:
                value = getattr(record, attribute)
                row.append(value)
            writer.writerow(row)

def main(json_dir,table_out):
    json_files = glob.glob(json_dir+"/*.json")
    records=[]
    for json_file in json_files:
        record = ProteinRecord.from_json(json_file)
        records.append(record)
    write_metadata_summary(records, table_out)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("json_dir", type=str, help="Path to directory with json files")
    parser.add_argument("table_out", type=str, help="Path to output table in .csv format")
    args = parser.parse_args()
    #Run the main script
    main(args.json_dir, args.table_out)