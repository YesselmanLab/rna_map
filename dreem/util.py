import os
import subprocess
import shutil

# TODO remove too much requirements
from Bio import SeqIO


def safe_rmdir(dir_name):
    if os.path.isdir(dir_name):
        shutil.rmtree(dir_name)


def safe_mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


def fasta_to_dict(fasta_file):
    """
    Parse a FASTA file
    Args:
        fasta_file (string): Path to FASTA file
    Returns:
        refs_seq (dict): Sequences of the ref genomes in the file
    """

    refs_seq = {}
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            refs_seq[record.id] = str(record.seq)
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
