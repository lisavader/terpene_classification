import argparse
import glob
import os
import re
from collections import defaultdict
from subprocess import Popen, PIPE
from io import StringIO
from Bio import SearchIO

from fasta_parsing import read_fasta

def main(hmm_database, fasta, accessions_path, scores_out):
    accessions_all = read_fasta(fasta).keys()
    hmms_accessions = defaultdict(lambda: defaultdict(dict))

    for accessions_file in glob.glob(accessions_path):
        hmm_name = os.path.splitext(os.path.basename(accessions_file))[0]
        with open(accessions_file,"r") as stream:
            accessions_in = set()
            for line in stream:
                accession = re.sub("C-|N-","",line.rstrip())
                assert accession in accessions_all, f"Accession {accession} in {accessions_file} is not present in fasta file {fasta}"
                accessions_in.add(accession)

        accessions_out = list(set(accessions_all) - accessions_in)
        accessions_in = list(accessions_in)
        hmms_accessions[hmm_name]["accessions_in"] = accessions_in
        hmms_accessions[hmm_name]["accessions_out"] = accessions_out

    hmmscan_results = run_hmmscan(hmm_database, fasta)

    all_results = defaultdict(dict)
    for result in hmmscan_results:
        for hsp in result.hsps:
            try:
                all_results[hsp.hit_id][hsp.query_id].append(hsp.bitscore)
            except:
                all_results[hsp.hit_id][hsp.query_id] = [hsp.bitscore]

    scores_per_hmm = defaultdict(lambda: defaultdict(dict))
    for hmm_name in hmms_accessions:
        for accession in hmms_accessions[hmm_name]["accessions_in"]:
            try:
                #Retrieve best score for hmm with the in_group sequence left out
                topscore = max(all_results[hmm_name+"_"+accession][accession])
            except:
                topscore = 0
            scores_per_hmm[hmm_name][accession]["score"] = topscore
            scores_per_hmm[hmm_name][accession]["ingroup"] = True

        for accession in hmms_accessions[hmm_name]["accessions_out"]:
            try:
                #Retrieve best score for the general hmm
                topscore = max(all_results[hmm_name][accession])
            except:
                topscore = 0
            scores_per_hmm[hmm_name][accession]["score"] = topscore
            scores_per_hmm[hmm_name][accession]["ingroup"] = False

    with open(scores_out, 'w') as stream:
        stream.write('\t'.join(["hmm_profile", "accession", "score", "ingroup"])+'\n')
        for hmm_name in scores_per_hmm:
            for accession in scores_per_hmm[hmm_name]:
                score = scores_per_hmm[hmm_name][accession]["score"]
                ingroup = scores_per_hmm[hmm_name][accession]["ingroup"]
                stream.write('\t'.join(map(str,[hmm_name, accession, score, ingroup]))+'\n')

def run_hmmscan(hmm_database, fasta):
    proc = Popen(["hmmscan", "--noali", "--domT", "0", hmm_database, fasta], stdout=PIPE)
    hmmscan_out = proc.communicate()[0].decode()
    results = SearchIO.parse(StringIO(hmmscan_out), 'hmmer3-text')
    return results

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("hmm_database", type=str, help="Database of hmm files to scan against")
    parser.add_argument("fasta", type=str, help="Fast file with (full-length) input sequences")
    parser.add_argument("accessions_dir", type=str, help="Directory with accessions per hmm")
    parser.add_argument("scores_out", type=str, help="Output file containing all bitscores")
    args = parser.parse_args()
    #Run the main script
    main(**vars(args))
