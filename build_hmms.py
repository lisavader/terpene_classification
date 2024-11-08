import argparse
import glob
import os
import logging
import re
from subprocess import Popen, PIPE

from fasta_parsing import read_fasta

def main(accessions_path, hmm_dir, fasta, leave_one_out, force):
    fasta_dict = read_fasta(fasta)
    for accessions_file in glob.glob(accessions_path):
        hmm_name = os.path.splitext(os.path.basename(accessions_file))[0]
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
        if len(sequences) != len(accessions):
            logging.warning(f"Skipping {hmm_name} and all LOO hmms; {str(len(accessions))} sequences are required but only {str(len(sequences))} found.")
            continue
        build_hmm(sequences, hmm_dir, hmm_name, force)

        #Build leave-one-out hmms
        if leave_one_out:
            base_accessions = list(set([re.sub("C-|N-","",accession) for accession in accessions]))
            for base_accession in base_accessions:
                sequences_loo = dict(sequences)
                for accession in sequences.keys():
                    if accession.endswith(base_accession):
                        del sequences_loo[accession]
                hmm_name_loo = hmm_name+"_"+base_accession
                build_hmm(sequences_loo, hmm_dir, hmm_name_loo, force)

def build_hmm(sequences, hmm_dir, hmm_name, force):
    hmm_file = os.path.join(hmm_dir, hmm_name+".hmm")
    if os.path.isfile(hmm_file) and not force:
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
    parser.add_argument("accessions_path", type=str, help="Path to text files with accessions with which to build a hmm. May contain wildcards.")
    parser.add_argument("hmm_dir", type=str, help="Directory for storing the hmm files")
    parser.add_argument("fasta", type=str, help="Fasta file that contains the sequences for each accession")
    parser.add_argument("-l", "--leave_one_out", action="store_true", help="Build additional hmms with one of the sequences left out.")
    parser.add_argument("-v", "--verbose", action="store_const", dest="loglevel", const=logging.INFO, help="Enable verbose mode")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    #Run the main script
    main(args.accessions_path, args.hmm_dir, args.fasta, args.leave_one_out, args.force)