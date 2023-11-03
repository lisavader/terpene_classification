import argparse
import glob
import ast

from Records import ProteinRecord

def write_annotation_labels(records,file_path):
    with open(file_path,"w") as file:
        file.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for record in records:
            file.write(record.accession+"\t"+record.protein_name+"\n")

def write_annotation_colours(records,file_path,group_by,group_colours):
    with open(group_colours, "r") as grouping_data:
        grouping_dict = ast.literal_eval(grouping_data.read())
    with open(file_path,"w") as file:
        file.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for record in records:
            colour = ""
            attribute_value = getattr(record, group_by)
            if isinstance(attribute_value, list):
                for item in attribute_value:
                    if item in grouping_dict:
                        colour = grouping_dict.get(item)
            else:
                colour = grouping_dict.get(attribute_value)
            file.write(record.accession+"\t"+"label_background"+"\t"+colour+"\n")

def main(json_dir,annotation_labels,annotation_colours,group_by,group_colours):
    json_files = glob.glob(json_dir+"/*.json")
    records=[]
    for json_file in json_files:
        record = ProteinRecord.from_json(json_file)
        records.append(record)
    if annotation_labels:
        write_annotation_labels(records,annotation_labels)
    if annotation_colours:
        write_annotation_colours(records,annotation_colours,group_by,group_colours)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("json_dir", type=str, help="Path to directory with json files")
    parser.add_argument("--annotation_labels", type=str, help="Path to annotation labels output")
    parser.add_argument("--annotation_colours", type=str, help="Path to annotation colours output")
    parser.add_argument("--group_by", type=str, help="Which attribute to use for grouping")
    parser.add_argument("--group_colours", type=str, help="Dictionary showing how to colour each group")
    args = parser.parse_args()
    #Run the main script
    main(args.json_dir,args.annotation_labels,args.annotation_colours,args.group_by,args.group_colours)