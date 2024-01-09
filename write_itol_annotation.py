import argparse
import json
import sys

from fasta_parsing import read_fasta
from Records import ProteinRecord

def write_labels(records,file_path):
    with open(file_path,"w") as file:
        file.write("LABELS\nSEPARATOR COMMA\nDATA\n")
        for record in records:
            file.write(record.accession+","+record.protein_name+"\n")

def write_annotation_groups(records,file_path,group_by,annotation_type,group_values):
    with open(group_values,"r") as file:
        legend = json.load(file)

    with open(file_path,"w") as file:
        if annotation_type == "colour":
            file.write("TREE_COLORS\nSEPARATOR COMMA\nDATA\n")
        elif annotation_type == "shape":
            shapes = list(legend.values())
            shapes_string = ",".join(shapes)
            file.write("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,shapes\nCOLOR,#000000\nFIELD_SHAPES,"+shapes_string+
                       "\nFIELD_LABELS,"+shapes_string+"\nSYMBOL_SPACING,-23\nDATA\n")

        for record in records:
            try:
                attribute = str(getattr(record, group_by))
            except AttributeError:
                print("Error: No attribute with name '"+group_by+"' found in metadata of record "+record.accession)
                sys.exit(1)

            value = legend.get(attribute)
            if annotation_type == "colour":
                file.write(record.accession+","+"label_background"+","+value+"\n")
            elif annotation_type == "shape":
                shape_values = ["-1"] * len(shapes)
                index = shapes.index(value)
                shape_values[index] = "0"
                shape_values_string = ",".join(shape_values)
                file.write(record.accession+","+shape_values_string+"\n")

def main(fasta_in,json_dir,labels,group_annotation,group_by,annotation_type,group_values):
    json_files = []
    records=[]
    fasta_dict = read_fasta(fasta_in)
    for header in fasta_dict.keys():
        json_file = json_dir+header+".json"
        json_files.append(json_file)
    for json_file in json_files:
        record = ProteinRecord.from_json(json_file)
        records.append(record)
    if labels:
        write_labels(records,labels)
    if group_annotation:
        write_annotation_groups(records,group_annotation,group_by,annotation_type,group_values)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to input fasta")
    parser.add_argument("json_dir", type=str, help="Path to directory with json files")
    parser.add_argument("--labels", type=str, help="Path to leaf labels output")
    parser.add_argument("--group_annotation", type=str, help="Path to group annotation output")
    parser.add_argument("--group_by", type=str, help="Which attribute to use for grouping")
    parser.add_argument("--annotation_type", type=str, help="The type of group annotation (choose 'colour' or 'shape')")
    parser.add_argument("--group_values", type=str, help="Supply a json file with a colour/shape for each group.")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in,args.json_dir,args.labels,args.group_annotation,args.group_by,args.annotation_type,args.group_values)