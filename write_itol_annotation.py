import argparse
import json
import sys

from fasta_parsing import read_fasta
from Records import ProteinRecord

def write_annotation_labels(records,file_path):
    with open(file_path,"w") as file:
        file.write("LABELS\nSEPARATOR COMMA\nDATA\n")
        for record in records:
            file.write(record.accession+","+record.protein_name+"\n")

def write_annotation_colours(records,file_path,group_by,legend_out):
    colour_list = ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"]
    legend = {}
    i = 0
    with open(file_path,"w") as file:
        file.write("TREE_COLORS\nSEPARATOR COMMA\nDATA\n")
        for record in records:
            try:
                attribute = str(getattr(record, group_by))
            except AttributeError:
                print("Error: No attribute with name '"+group_by+"' found in metadata of record "+record.accession)
                sys.exit(1)

            if attribute not in legend.keys():
                colour = colour_list[i]
                legend[attribute] = colour
                i = i + 1
            else:
                colour = legend.get(attribute)
            file.write(record.accession+","+"label_background"+","+colour+"\n")
    if legend_out:
        with open(legend_out,"w") as file:
            json.dump(legend, file)

def main(fasta_in,json_dir,annotation_labels,annotation_colours,group_by,legend_out):
    json_files = []
    records=[]
    fasta_dict = read_fasta(fasta_in)
    for header in fasta_dict.keys():
        json_file = json_dir+header+".json"
        json_files.append(json_file)
    for json_file in json_files:
        record = ProteinRecord.from_json(json_file)
        records.append(record)
    write_annotation_labels(records,annotation_labels)
    if annotation_colours:
        write_annotation_colours(records,annotation_colours,group_by,legend_out)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_in", type=str, help="Path to input fasta")
    parser.add_argument("json_dir", type=str, help="Path to directory with json files")
    parser.add_argument("annotation_labels", type=str, help="Path to annotation labels output")
    parser.add_argument("--annotation_colours", type=str, help="Path to annotation colours output")
    parser.add_argument("--group_by", type=str, help="Which attribute to use for grouping")
    parser.add_argument("--legend_out", type=str, help="Path to legend output")
    args = parser.parse_args()
    #Run the main script
    main(args.fasta_in,args.json_dir,args.annotation_labels,args.annotation_colours,args.group_by,args.legend_out)