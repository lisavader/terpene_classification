import argparse

from InterproRecord import InterproRecord
from fasta_parsing import read_fasta

def curate_record(record):
    if record.protein_name == "Aristolochene_synthase_in_complex_with_12,13_Difluorofarnesyl_diphosphate":
        record.enzyme_class = ["terpene synthase"]
        record.enzyme_subclass = ["sesquiterpene synthase"]
    if "gibberellin_cluster-ggpp-synthase" in record.protein_name.lower() or record.protein_name in \
        ["Geranylgeranyl_pyrophosphate_synthetase,_putative_[includes:_dimethylallyltranstransferase_(Ec_2.5.1.1)_geranyltranstransferase_(Ec_2.5.1.10)_farnesyltranstransferas_(Ec_2.5.1.29)]",
        "Geranylgeranyl_pyrophosphate_synthase_(Fusicoccadiene_synthase)"]:
        record.enzyme_class = ["prenyltransferase"]
        record.enzyme_subclass = ["GGPP synthase"]
    return record

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
        record = curate_record(record)
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