import argparse
import urllib
import json
from bs4 import BeautifulSoup
from Bio import Entrez, GenBank, SeqIO
from io import StringIO

from Records import ProteinRecord

Entrez.email = "lisvad@dtu.dk"

def efetch_protein(accession, output_type):
    handle = Entrez.efetch(db = "protein", id = accession, rettype = output_type)
    result = handle.read()
    return result

def elink_nucleotide_to_protein(accession):
    handle = Entrez.elink(dbfrom = "nucleotide", db = "protein", id = accession)
    result = handle.read()
    return result

def download_genbank(accession):
    #search for accession in protein db
    try:
        genbank_string = efetch_protein(accession, "gb")
    #if not a protein accession, search in nucleotide db
    except urllib.error.HTTPError:
        elink_result = elink_nucleotide_to_protein(accession)
        soup = BeautifulSoup(elink_result, features="xml")
        if len(soup.LinkSetDb.find_all("Link")) > 1:
            raise ValueError("Accession: "+accession+" is linked to multiple protein ids")
        else:
            uid = soup.Link.Id.string
            genbank_string = efetch_protein(uid, "gb")

    return genbank_string

def write_metadata(seqrecord, shared_metadata, json_dir):
    database = "ncbi"
    header = seqrecord.id+" "+seqrecord.description
    proteinrecord = ProteinRecord(database=database, header=header, accession=seqrecord.id, reviewed=None,
                                  protein_name=seqrecord.description, taxid=None)
    #add organism information
    proteinrecord.organism_name = seqrecord.annotations["organism"]
    for division in ["Ascomycota", "Basidiomycota"]:
        if division in seqrecord.annotations["taxonomy"]:
            proteinrecord.organism_division = division
        else:
            proteinrecord.organism_division = "Other"
    #add shared metadata
    if shared_metadata:
        with open(shared_metadata,"r") as file:
            metadata_dict = json.load(file)
        for attribute, value in metadata_dict.items():
            setattr(proteinrecord, attribute, value)
    #assign enzyme type
    if not proteinrecord.enzyme_type:
        proteinrecord.assign_enzyme_type()
    #store ProteinRecord object as json
    file_path = json_dir+"/"+proteinrecord.accession+".json"
    json_object = proteinrecord.to_json()
    with open(file_path,"w") as file:
        file.write(json_object)

def main(accessions, fasta_out, json_dir, shared_metadata):
    with open(accessions,"r") as file:
        accessions = []
        for line in file:
            accessions.append(line.rstrip())
    seqrecords = []
    for accession in accessions:
        genbank_string = download_genbank(accession)
        seqrecord = SeqIO.read(StringIO(genbank_string), format="gb")
        seqrecords.append(seqrecord)
        write_metadata(seqrecord, shared_metadata, json_dir)
    with open(fasta_out,"w") as output_file:
        SeqIO.write(seqrecords, output_file, "fasta")

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", type=str, help="File with ncbi protein or nucleotide accessions, separated by new line")
    parser.add_argument("fasta_out", type=str, help="File to output fasta files")
    parser.add_argument("json_dir", type=str, help="Path to json directory for storing metadata")
    parser.add_argument("shared_metadata", type=str, nargs='?', default=None, help="Optional: Metadata that is shared between all sequences, in json format.")
    args = parser.parse_args()
    #Run the main script
    main(args.accessions, args.fasta_out, args.json_dir, args.shared_metadata)