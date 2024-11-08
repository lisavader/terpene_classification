import argparse
import json
import os

from fasta_parsing import read_fasta
from Records import ProteinRecord

def write_labels(json_files,file_path):
    with open(file_path,"w") as file:
        file.write("LABELS\nSEPARATOR COMMA\nDATA\n")
        for accession, json_file in json_files.items():
            record = ProteinRecord.from_json(json_file)
            file.write(accession+","+record.protein_name+"\n")

def write_annotation_groups(json_files,file_path,group_by,annotation_type,group_values,shapes_offset):
    with open(group_values,"r") as file:
        legend = json.load(file)

    with open(file_path,"w") as file:
        if annotation_type == "colour":
            file.write("TREE_COLORS\nSEPARATOR COMMA\nDATA\n")
        elif annotation_type == "colour_strip":
            file.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,colour_strip\nCOLOR,#000000\nDATA\n")
        elif annotation_type == "shape":
            shapes = list(legend.values())
            shapes_string = ",".join(shapes)
            file.write("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,shapes\nCOLOR,#000000\nFIELD_SHAPES,"+shapes_string+
                       "\nFIELD_LABELS,"+shapes_string+"\nSYMBOL_SPACING, "+str(shapes_offset)+"\nDATA\n")

        for accession, json_file in json_files.items():
            record = ProteinRecord.from_json(json_file)
            try:
                attribute = getattr(record, group_by)
            except AttributeError:
                print("Error: No attribute with name '"+group_by+"' found in metadata of record "+accession)
                continue

            if type(attribute) == list:
                attribute = str(attribute)
            try:
                value = legend[attribute]
            except KeyError:
                print("Error: Attribute \'"+str(attribute)+"\' of record "+accession+" not found in group_values file.")
                continue

            if annotation_type == "colour":
                file.write(accession+","+"label_background"+","+value+"\n")
            elif annotation_type == "colour_strip":
                file.write(accession+","+value+"\n")
            elif annotation_type == "shape":
                shape_values = ["-1"] * len(shapes)
                index = shapes.index(value)
                shape_values[index] = "0"
                shape_values_string = ",".join(shape_values)
                file.write(accession+","+shape_values_string+"\n")

def main(fasta_in,json_dir,labels,group_annotation,group_by,annotation_type,group_values,shapes_offset):
    json_files = {}
    fasta_dict = read_fasta(fasta_in)
    for header in fasta_dict.keys():
        if os.path.isfile(json_dir+header+".json"):
            json_file = json_dir+header+".json"
        else:
            json_file = json_dir+header[2:]+".json"
        accession = header
        json_files[accession] = json_file
    if labels:
        write_labels(json_files,labels)
    if group_annotation:
        write_annotation_groups(json_files,group_annotation,group_by,annotation_type,group_values,shapes_offset)

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
    parser.add_argument("--shapes_offset", type=int, help="Distance between the different shapes (an itol artifact)."+
                        "Choose a suitable negative offset to make them appear on the same line.")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in,args.json_dir,args.labels,args.group_annotation,args.group_by,args.annotation_type,args.group_values,args.shapes_offset)