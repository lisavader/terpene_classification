import argparse
import urllib
from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
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

def write_metadata(seqrecord, enzyme_type, json_dir):
    database = "ncbi"
    header = seqrecord.id+" "+seqrecord.description
    proteinrecord = ProteinRecord(database=database, header=header, accession=seqrecord.id, reviewed=None,
                                  protein_name=seqrecord.description, taxid=None)
    #add organism information
    proteinrecord.organism_name = seqrecord.annotations["organism"]
    for clade in ["Ascomycota", "Basidiomycota", "Bacteria", "Viridiplantae"]:
        if clade in seqrecord.annotations["taxonomy"]:
            proteinrecord.organism_category = clade
            break
        proteinrecord.organism_category = "Other"
    #add enzyme type
    if not enzyme_type == "none":
        proteinrecord.enzyme_type = enzyme_type.rstrip().split(", ")
    else:
        proteinrecord.enzyme_type = None
    #all manually added sequences are reviewed
    proteinrecord.reviewed = True
    #store ProteinRecord object as json
    file_path = json_dir+"/"+proteinrecord.accession+".json"
    json_object = proteinrecord.to_json()
    with open(file_path,"w") as file:
        file.write(json_object)

def main(ncbi_file, fasta_out, json_dir):
    accessions = {}
    seqrecords = []
    with open(ncbi_file,"r") as file:
        next(file)  #skip header
        for line in file:
            accession, enzyme_type = line.split('\t')[0:2]
            accessions[accession] = enzyme_type
    for accession, enzyme_type in accessions.items():
        genbank_string = download_genbank(accession)
        seqrecord = SeqIO.read(StringIO(genbank_string), format="gb")
        write_metadata(seqrecord, enzyme_type, json_dir)
        seqrecords.append(seqrecord)
    with open(fasta_out,"w") as output_file:
        SeqIO.write(seqrecords, output_file, "fasta")

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("ncbi_file", type=str, help="Tsv file with ncbi protein or nucleotide accessions in first column, enzyme type in second column")
    parser.add_argument("fasta_out", type=str, help="File to output fasta files")
    parser.add_argument("json_dir", type=str, help="Path to json directory for storing metadata")
    args = parser.parse_args()
    #Run the main script
    main(**vars(args))