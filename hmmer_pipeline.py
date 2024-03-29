import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from run_hmmscan import run_hmmscan

def queryresult_to_df(queryresult):
    #inspired by: https://stackoverflow.com/questions/62012615/convert-a-hmmer-tblout-output-to-a-pandas-dataframe
    hit_attributes = ["id","evalue","bitscore","bias"]
    hits = defaultdict(list)
    for result in queryresult:
        for hit in result.hits:
            query_id = getattr(result, "id").split('|',1)[0]
            hits["query_id"].append(query_id)
            for attribute in hit_attributes:
                hits[attribute].append(getattr(hit, attribute))
    hmmer_df = pd.DataFrame.from_dict(hits)
    return hmmer_df

def main(hmmfile, fasta, hmmdetails, output, correct_accessions, opts):
    if opts:
        opts_list = opts.split(',')
    else:
        opts_list = None
    queryresult = run_hmmscan(hmmfile, fasta, opts_list)
    hmmer_df = queryresult_to_df(queryresult)
    if hmmdetails:
        details = pd.read_csv(hmmdetails, delimiter='\t', header=None, index_col=0)
        cutoffs = details.iloc[:,1]
        cutoff_dict = cutoffs.to_dict()
        hmmer_df["cutoff_score"] = hmmer_df["id"].apply(lambda x: cutoff_dict.get(x))
        hmmer_df["significant"] = np.where(hmmer_df["bitscore"] >= hmmer_df["cutoff_score"], "True", "False")
    if correct_accessions:
        with open(correct_accessions,"r") as file:
            accessions = []
            for line in file:
                accessions.append(line.rstrip())
        hmmer_df["correct"] = hmmer_df["query_id"].isin(accessions)
    hmmer_df.to_csv(output, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hmmfile", type=str, help="the path to a HMM file to use in scanning")
    parser.add_argument("fasta", type=str, help="a string containing input sequences in fasta format")
    parser.add_argument("--hmmdetails", type=str, help="Optional: Path to the antiSMASH file hmmdetails.txt describing cutoffs for each profile."+
                        "This is used to calculate if hits are 'significant' according to antiSMASH (above the cut-off)")
    parser.add_argument("output", type=str, help="Path to output .tsv file")
    parser.add_argument("--correct_accessions", type=str, help="Optional: supply a file with correct accessions that are supposed"+
                        " to be picked up by the hmmm profile(s)")
    parser.add_argument("--opts", type=str, help="Add options to hmmscan (as string where fields are separated by commas)")
    args = parser.parse_args()
    main(args.hmmfile, args.fasta, args.hmmdetails, args.output, args.correct_accessions, args.opts)