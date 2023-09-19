import pandas as pd
import numpy as np
import argparse
import json
from Bio import SearchIO
from collections import defaultdict

def parse_hmmer_tbl(hmmer_tbl):
    #inspired by: https://stackoverflow.com/questions/62012615/convert-a-hmmer-tblout-output-to-a-pandas-dataframe
    attributes = ["id","evalue","bitscore","bias"]
    hits = defaultdict(list)
    with open(hmmer_tbl) as handle:
        for result in SearchIO.parse(handle, "hmmer3-tab"):
            for hit in result.hits:
                hits["query_id"].append(getattr(result, "id"))
                for attribute in attributes:
                    hits[attribute].append(getattr(hit, attribute))
    hmmer_df = pd.DataFrame.from_dict(hits)
    return hmmer_df

def main(hmmer_tbl,hmm_details,json_dir,output):
    hmmer_df = parse_hmmer_tbl(hmmer_tbl)
    details = pd.read_csv(hmm_details, delimiter='\t', header=None, index_col=0)
    cutoffs = details.iloc[:,1]
    cutoff_dict = cutoffs.to_dict()
    hmmer_df["cutoff_score"] = hmmer_df["id"].apply(lambda x: cutoff_dict.get(x))
    hmmer_df["significant"] = np.where(hmmer_df["bitscore"] >= hmmer_df["cutoff_score"], "True", "False")
    for attribute in ["protein_name","enzyme_class","enzyme_subclass","organism_name","organism_division","review_status"]:
        hmmer_df[attribute] = hmmer_df["query_id"].apply(lambda x: json.load(open(json_dir+"/"+x+".json")).get(attribute))
    hmmer_df.to_csv(output, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hmmer_tbl", type=str, help="Path to a HMMer output file in tabular format")
    parser.add_argument("hmm_details", type=str, help="Path to the antiSMASH file hmmdetails.txt describing cutoffs for each profile")
    parser.add_argument("json_dir", type=str, help="Path to directory with json files containing extra info for each query sequence")
    parser.add_argument("output", type=str, help="Path to output .tsv file")
    args = parser.parse_args()
    main(args.hmmer_tbl, args.hmm_details, args.json_dir, args.output)