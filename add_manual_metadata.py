import argparse
import json
from Records import ProteinRecord

def update_json(metadata, json_dir):
    file_path = json_dir+metadata["accession"]+".json"
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        print("File "+file_path+" not found!\nCreating .json file at "+file_path)
        data = {}
        data["accession"] = metadata["accession"]
    for key, value in list(metadata.items())[1:]:
        if value:
            data[key] = value
        else:
            data[key] = "none"
    with open(file_path,"w") as file:
        json.dump(data, file, default=lambda o: o.__dict__, indent=4)

def main(data_in, json_dir):
    with open(data_in, "r", encoding='utf-16') as csv_file:
        header = next(csv_file)
        column_names = header.strip('\n').split('\t')
        for line in csv_file:
            line = line.strip('\n')
            fields = line.split('\t')
            metadata = dict(zip(column_names,fields))
            metadata["enzyme_type"] = metadata["enzyme_type"].split(',')    #convert to list
            update_json(metadata, json_dir)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("data_in", type=str, help="Path to .tsv file containing metadata")
    parser.add_argument("json_dir", type=str, help="Path to json output directory")
    args = parser.parse_args()
    #Run the main script
    main(args.data_in, args.json_dir)