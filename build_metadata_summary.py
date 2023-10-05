import argparse
import glob
import csv
import json

def main(json_dir, output):
    json_files = glob.glob(json_dir+"/*.json")
    with open(output, 'w') as csv_out:
        writer = csv.writer(csv_out)
        header = ["accession","protein_name","enzyme_class","enzyme_subclass","organism_name","organism_division","review_status"]
        writer.writerow(header)
        for json_file in json_files:
            row = []
            for attribute in header:
                value = json.load(open(json_file)).get(attribute)
                row.append(value)
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("json_dir", type=str, help="Path to directory with json files containing extra info for each query sequence")
    parser.add_argument("output", type=str, help="Path to output .csv file")
    args = parser.parse_args()
    main(args.json_dir, args.output)