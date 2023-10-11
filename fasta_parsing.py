import logging
from typing import Dict, List

# These functions are (slightly) modified from AntiSMASH v.7.0.1 (antismash/antismash/common/fasta.py)

def write_fasta(names: List[str], seqs: List[str], filename: str) -> None:
    """ Writes name/sequence pairs to file in FASTA format

        Argumnets:
            names: a list of sequence identifiers
            seqs: a list of sequences as strings
            filename: the filename to write the FASTA formatted data to

        Returns:
            None
    """
    with open(filename, "w", encoding="utf-8") as out_file:
        for name, seq in zip(names, seqs):
            out_file.write(f">{name}\n{seq}\n")


def read_fasta(filename: str) -> Dict[str, str]:
    """ Reads a fasta file into a dictionary

        Arguments:
            filename: the path to the FASTA file to read

        Returns:
            a dictionary mapping sequence ID to sequence

    """
    ids = []
    sequence_info = []
    with open(filename, "r", encoding="utf-8") as fasta:
        current_seq: List[str] = []
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                ids.append(line[1:].replace(" ", "_"))
                if current_seq:
                    sequence_info.append("".join(current_seq))
                    current_seq.clear()
            else:
                if not ids:
                    raise ValueError("Sequence before identifier in fasta file")
                if not line.replace("-", "z").isalpha():
                    raise ValueError("Sequence contains non-alphabetic characters")
                current_seq.append(line)
    if current_seq:
        sequence_info.append("".join(current_seq))
    if len(ids) != len(sequence_info):
        raise ValueError("Fasta files contains different counts of sequences and ids")
    if not ids:
        logging.debug("Fasta file %s contains no sequences", filename)
        raise ValueError("Fasta file contains no sequences")
    return dict(zip(ids, sequence_info))