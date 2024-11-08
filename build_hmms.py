import argparse
import glob
import os
import logging
import re
from subprocess import Popen, PIPE

from fasta_parsing import read_fasta

def main(accessions_dir, hmm_dir, fasta, leave_one_out):
    fasta_dict = read_fasta(fasta)
    for accessions_file in glob.glob(os.path.join(accessions_dir,"*")):
        hmm_name = os.path.basename(accessions_file).strip(".txt")
        with open(accessions_file,"r") as stream:
            accessions = []
            for line in stream:
                accessions.append(line.rstrip())

        #Build main hmm
        sequences = {}
        for accession in accessions:
            for header,seq in fasta_dict.items():
                if header == accession:
                    sequences[header] = seq
        build_hmm(sequences, hmm_dir, hmm_name)

        #Build leave-one-out hmms
        if leave_one_out:
            base_accessions = list(set([re.sub("C-|N-","",accession) for accession in accessions]))
            for base_accession in base_accessions:
                sequences_loo = dict(sequences)
                for accession in sequences.keys():
                    if accession.endswith(base_accession):
                        del sequences_loo[accession]
                hmm_name_loo = hmm_name+"_"+base_accession
                build_hmm(sequences_loo, hmm_dir, hmm_name_loo)

def build_hmm(sequences, hmm_dir, hmm_name):
    hmm_file = os.path.join(hmm_dir, hmm_name+".hmm")
    if os.path.isfile(hmm_file):
        logging.info(f"Skipping {hmm_name}, file already exists at {hmm_file}")
        return
    else:
        sequence_string = ""
        for header,seq in sequences.items():
            sequence_string += f">{header}\n{seq}\n"
        alignment = build_alignment(sequence_string)
        logging.info(f"Writing hmm {hmm_name} based on {str(len(sequences))} sequences to {hmm_file}")
        run_hmmbuild(alignment, hmm_name, hmm_file)
    return

def build_alignment(sequence_string):
    proc = Popen(["muscle", "-quiet"], stdin=PIPE, stdout=PIPE)
    alignment = proc.communicate(input=sequence_string.encode())[0].decode()
    return alignment

def run_hmmbuild(alignment, hmm_name, hmm_file):
    proc = Popen(["hmmbuild", "-n", hmm_name, "--informat", "afa", hmm_file, "-"], stdin=PIPE, stdout=PIPE)
    proc.communicate(input=alignment.encode())
    return

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions_dir", type=str, help="Directory that contains text files with accessions with which to build a hmm.")
    parser.add_argument("hmm_dir", type=str, help="Directory for storing the hmm files")
    parser.add_argument("fasta", type=str, help="Fasta file that contains the sequences for each accession")
    parser.add_argument("--leave_one_out", action="store_true", help="Build additional hmms with one of the sequences left out.")
    parser.add_argument("-v", "--verbose", action="store_const", dest="loglevel", const=logging.INFO, help="Enable verbose mode")
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    #Run the main script
    main(args.accessions_dir, args.hmm_dir, args.fasta, args.leave_one_out)