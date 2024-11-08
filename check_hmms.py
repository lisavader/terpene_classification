import argparse
import glob
import os
import statistics
from collections import defaultdict
from subprocess import Popen, PIPE
from io import StringIO
from Bio import SearchIO

from fasta_parsing import read_fasta

def main(hmm_database, fasta, accessions_dir, scores_out, stats_out):
    accessions_all = read_fasta(fasta).keys()
    hmmscan_results = run_hmmscan(hmm_database, fasta)

    all_results = defaultdict(dict)
    for result in hmmscan_results:
        for hsp in result.hsps:
            try:
                all_results[hsp.hit_id][hsp.query_id].append(hsp.bitscore)
            except:
                all_results[hsp.hit_id][hsp.query_id] = [hsp.bitscore]

    grouped_results = defaultdict(dict)
    for accessions_file in glob.glob(os.path.join(accessions_dir,"*")):
        hmm_name = os.path.basename(accessions_file).strip(".txt")
        with open(accessions_file,"r") as stream:
            accessions_in = []
            for line in stream:
                accessions_in.append(line.rstrip())
        accessions_out = list(set(accessions_all) - set(accessions_in))

        grouped_results[hmm_name] = {
            "ingroup": [],
            "outgroup": []
        }
        for accession in accessions_in:
            try:
                #Retrieve best score for hmm with the in_group sequence left out
                topscore = max(all_results[hmm_name+"_"+accession][accession])
            except:
                topscore = 0
            grouped_results[hmm_name]["ingroup"].append(topscore)
        for accession in accessions_out:
            try:
                #Retrieve best score for the general hmm
                topscore = max(all_results[hmm_name][accession])
            except:
                topscore = 0
            grouped_results[hmm_name]["outgroup"].append(topscore)

    with open(scores_out, 'w') as stream:
        stream.write('\t'.join(["hmm_profile", "sequence", "bitscores"])+'\n')
        for hmm_profile in all_results:
            for query, bitscores in all_results[hmm_profile].items():
                stream.write('\t'.join([hmm_profile, query, ','.join(map(str, bitscores))])+'\n')

    with open(stats_out, 'w') as stream:
        stream.write('\t'.join(["hmm", "median_in", "min_in", "max_in", "median_out", "min_out", "max_out", "proposed_cutoff"])+'\n')
        for hmm_name in grouped_results:
            median_in = statistics.median(grouped_results[hmm_name]["ingroup"])
            min_in = min(grouped_results[hmm_name]["ingroup"])
            max_in = max(grouped_results[hmm_name]["ingroup"])
            median_out = statistics.median(grouped_results[hmm_name]["outgroup"])
            min_out = min(grouped_results[hmm_name]["outgroup"])
            max_out = max(grouped_results[hmm_name]["outgroup"])
            proposed_cutoff = int(statistics.mean([min_in, max_out]))
            stream.write('\t'.join(map(str,[hmm_name, median_in, min_in, max_in, median_out, min_out, max_out, proposed_cutoff]))+'\n')

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
    parser.add_argument("stats_out", type=str, help="Output file with bitscore statistics per hmm profile")
    args = parser.parse_args()
    #Run the main script
    main(**vars(args))
