import argparse
import re
import ast
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
from ete3 import NCBITaxa

#connector to ncbi database
ncbi = NCBITaxa()

#build query lists
monoterpene=["monoterpene","iridoid"]
sesquiterpene=["sesquiterpene","sesquiterpenoid","trichodiene","trichodience","protoilludene","aristolochene","guaiadiene","germacrene","farnesene","viridiflorene","cuprenene"
               ,"muurolene","longiborneol","presilphiperfolanol","eremophilene","cadinene","santalene","isozizaene","pentalenene","pristinol"]
diterpene=["diterpene","variediene","fusicoccadiene","conidiogenone","araneosene","kaurene","phomopsene","dolastadiene","aphidicolanol","paxilline","BcPAX1","copalyl diphosphate"
            ,"gibberellin","phyllocladanol"]
sesterterpene=["sesterterpene","ophiobolin","preasperterpenoid","preaspterpenacid","stellatatriene","astellifadiene","quiannulatene","sesterbrasiliatriene","sesterfisherol"
                ,"preterpestacin","aspergildiene","fusoxypene"]
triterpene=["triterpene","macrophomene","talaropentaene"]
tetraterpene=["tetraterpene","lycopene"]

FPP=["FPP","farnesyl pyrophosphate","farnesyl diphosphate","dimethylallyltranstransferase","ERG20","geranyltranstransferase"]
GGPP=["GGPP","geranyl geranyl pyrophosphate","geranyl geranyl diphosphate","BTS1","farnesyltranstransferase"]
phytoene=["phytoene"]
squalene=["squalene","ERG9"]
HPP=["HPP","hexaprenyl pyrophosphate","hexaprenyl diphosphate","COQ1"]

generic_prenyltransferase=["prenyl transferase","polyprenyl","polyprenyl diphosphate"]
generic_terpene=["terpene","terpenoid"]

#convert to regular expressions
for query_list in [monoterpene,sesquiterpene,diterpene,sesterterpene,triterpene,tetraterpene,FPP,GGPP,HPP,generic_prenyltransferase]:
    for i,query in enumerate(query_list):
        regex = query
        if i != 0:      #skip for first query (is not a compound name)
            #compound suffixes -diene, -ene and -ol can be preceded by structural denotions between hyphens, e.g. dolasta-1(15),8-diene
            for suffix in ["diene","ene","ol"]:
                if suffix in query:
                    regex=re.sub(suffix+"$","(-.*-)?"+suffix,query)
                    break
        #spaces between words are optional, and sometimes replaced by a hyphen
        regex=re.sub(" ","[ -]?",regex)
        #replace original query with regex query
        query_list[i] = regex

class InterproRecord(SeqRecord):
    def __init__(self, seq_record):
        super().__init__(seq_record.seq, id=seq_record.id, description=seq_record.description)

        #split header into relevant fields
        self.header = self.description
        self.accession = self.header.split("|")[0]
        self.review_status = self.header.split("|")[1]
        self.protein_name = self.header.split("|")[2]
        self.taxid = self.header.split("|")[3].strip("taxID:")

        #rename id and description
        self.id = self.accession
        self.description = self.accession

        #initiate labels
        self.organism_name = ""
        self.organism_division = ""
        self.enzyme_class = []      #can have multiple classes and subclasses
        self.enzyme_subclass = []

    def assign_taxonomic_labels(self):
        taxid_translator = ncbi.get_taxid_translator([self.taxid])
        self.organism_name = list(taxid_translator.values())[0]     #the first taxonomic name found for this taxid
        lineage = ncbi.get_lineage(self.taxid)
        if 5204 in lineage:
            self.organism_division = "Basidiomycota"
        elif 4890 in lineage:
            self.organism_division = "Ascomycota"
        else:
            self.organism_division = "Other"

    def assign_enzyme_labels(self):
        # Recognise terpene synthases
        for query_list in [monoterpene,sesquiterpene,diterpene,sesterterpene,triterpene,tetraterpene]:
            for query in query_list:
                match = re.search(query,self.protein_name, flags=re.IGNORECASE)
                if match:
                    self.enzyme_class.append("terpene synthase")
                    self.enzyme_subclass.append(query_list[0]+" synthase")

        for query in generic_terpene:
            match = re.search(query,self.protein_name, flags=re.IGNORECASE)
            if match:
                self.enzyme_class.append("terpene synthase")

        # Recognise prenyltransferases
        for query_list in [FPP,GGPP,phytoene,squalene,HPP]:
            for query in query_list:
                match = re.search(query,self.protein_name, flags=re.IGNORECASE)
                if match:
                    self.enzyme_class.append("prenyltransferase")
                    self.enzyme_subclass.append(query_list[0]+" synthase")

        for query in generic_prenyltransferase:
            match = re.search(query,self.protein_name, flags=re.IGNORECASE)
            if match:
                self.enzyme_class.append("prenyltransferase")

        if not self.enzyme_class:
            self.enzyme_class.append("unknown")

        if not self.enzyme_subclass:
            self.enzyme_subclass.append("unknown")

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, indent=4)

def deduplicate_records(records):
    seq_set=set()
    deduplicated_records=[]

    for record in records:
        seq = record.seq
        if seq in seq_set:
            continue
        else:
            seq_set.add(seq)
            deduplicated_records.append(record)

    return deduplicated_records

def summarise_records(records, attribute_name, print_summary=True):
    #collect attributes values
    value_list = []
    for record in records:
        attribute_value = getattr(record, attribute_name)
        if attribute_value is set() or "":
            attribute_value = "None"
        value_list.append(attribute_value)
    #count occurences
    counter = Counter(value_list)
    ordered_counter = counter.most_common()

    if print_summary == True:
        #print line by line
        for item in ordered_counter:
            print(item[0]+' : '+str(item[1]))

    return ordered_counter

def select_records(records,selection_dict):
    selected_records=[]

    for record in records:
        for attribute_name,selection_values in selection_dict.items():
            match = False   #whether a match is found for this attribute
            attribute_value = getattr(record, attribute_name)
            if isinstance(attribute_value, list):    #check for all list items
                for item in attribute_value:
                    if item in selection_values:
                        match = True
                        break
            else:
                if attribute_value in selection_values:
                    match = True
            if match == False:  #all attributes in the selection dict must yield a match
                break
        if match == True:
            selected_records.append(record)
    return selected_records

def write_fasta(records,file_path):
    with open(file_path,"w") as file:
        SeqIO.write(records,file, "fasta")

def write_annotation_labels(records,file_path):
    with open(file_path,"w") as file:
        file.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for record in records:
            file.write(record.id+"\t"+record.protein_name+"\n")

def write_annotation_colours(records,file_path,group_by,group_colours):
    with open(group_colours, "r") as grouping_data:
        grouping_dict = ast.literal_eval(grouping_data.read())
    with open(file_path,"w") as file:
        file.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for record in records:
            colour = None
            attribute_value = getattr(record, group_by)
            if isinstance(attribute_value, set):
                for item in attribute_value:
                    if item in grouping_dict:
                        colour = grouping_dict.get(item)
            else:
                colour = grouping_dict.get(attribute_value)
            file.write(record.id+"\t"+"label_background"+"\t"+colour+"\n")

def write_json(records, json_dir):
    for record in records:
        file_path = json_dir+"/"+record.accession+".json"
        json_object = record.to_json()
        with open(file_path,"w") as file:
            file.write(json_object)

def main(input,output,selection,annotation_labels,annotation_colours,group_by,group_colours,json_dir):
    all_records = []
    for seq_record in SeqIO.parse(input, "fasta"):
        protein_record = InterproRecord(seq_record)
        protein_record.assign_enzyme_labels()
        protein_record.assign_taxonomic_labels()
        all_records.append(protein_record)
    with open(selection, "r") as selection_data:
        selection_dict = ast.literal_eval(selection_data.read())
    selected_records=select_records(all_records,selection_dict)
    write_fasta(selected_records, output)
    if annotation_labels:
        write_annotation_labels(selected_records,annotation_labels)
    if annotation_colours:
        write_annotation_colours(selected_records,annotation_colours,group_by,group_colours)
    if json_dir:
        write_json(all_records,json_dir)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="Path to fasta file with sequences from interpro")
    parser.add_argument("output", type=str, help="Path to output fasta file")
    parser.add_argument("selection", type=str, help="Path to file with selection criteria in dict format")
    parser.add_argument("--annotation_labels", type=str, help="Path to annotation labels output")
    parser.add_argument("--annotation_colours", type=str, help="Path to annotation colours output")
    parser.add_argument("--group_by", type=str, help="Which attribute to use for grouping")
    parser.add_argument("--group_colours", type=str, help="Legend showing how to colour each group")
    parser.add_argument("--json_dir", type=str, help="Path to json output dir")
    args = parser.parse_args()
    #Run the main script
    main(args.input,args.output,args.selection,args.annotation_labels,args.annotation_colours,args.group_by,args.group_colours,args.json_dir)