import os
import shutil
from pathlib import Path


def fasta_to_dict(fasta_file):
    """
    Parse a FASTA file
    :param fasta_file: (string): Path to FASTA file
    :return: refs_seq (dict): Sequences of the ref genomes in the file
    """

    refs_seq = {}
    with open(fasta_file, "r") as handle:
        for i, line in enumerate(handle.readlines()):
            if line[0] == ">":
                ref_name = line[1:].strip()
                refs_seq[ref_name] = ""
            else:
                refs_seq[ref_name] += line.strip()

    return refs_seq


def parse_phred_qscore_file(qscore_filename):
    phred_qscore = {}
    qscore_file = open(qscore_filename)
    qscore_file.readline()  # Ignore header line
    for line in qscore_file:
        line = line.strip().split()
        score, symbol = int(line[0]), line[1]
        phred_qscore[symbol] = score
    qscore_file.close()
    return phred_qscore


def get_filename(path):
    """
    get the filename from a path
    :param path: the path
    :return: the filename
    """
    path = Path(Path(path).stem)
    return str(path).rstrip(''.join(path.suffixes))