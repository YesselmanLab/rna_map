"""
run main program of rna_map
"""

import re
import os
import pandas as pd
import yaml

from rna_map.bit_vector import BitVectorGenerator

from rna_map.exception import DREEMInputException
from rna_map.logger import get_logger
from rna_map.parameters import (
    Inputs,
    get_default_params,
    validate_parameters,
)
from rna_map.mapping import Mapper
from rna_map.util import fasta_to_dict

log = get_logger("RUN")


# validate inputs #############################################################


def validate_fasta_file(fa: str) -> None:
    """
    ensure that the fasta file is in the correct format
    :param fa: path to the fasta file
    """
    with open(fa, "r", encoding="utf8") as f:
        lines = f.readlines()
    if len(lines) > 2000:
        log.warning(
            "fasta file contains more than 1000 sequences, this may lead"
            " to file generation issues. Its recommended to use --summary-output-only "
        )
    num = 0
    for i, l in enumerate(lines):
        l = l.rstrip()
        if len(l) == 0:
            raise DREEMInputException(
                f"blank line found on ln: {i}. These are not allowed in fastas."
            )
        # should be a reference sequence declartion
        if i % 2 == 0:
            num += 1
            if not l.startswith(">"):
                raise DREEMInputException(
                    f"reference sequence names are on line zero and even numbers."
                    f" line {i} has value which is not correct format in the fasta"
                )
            if l.startswith("> "):
                raise DREEMInputException(
                    f"there should be no spaces between > and reference name."
                    f"this occured on ln: {i} in the fasta file"
                )
        elif i % 2 == 1:
            if l.startswith(">"):
                raise DREEMInputException(
                    f"sequences should be on are on odd line numbers."
                    f" line {i} has value which is not correct format in fasta file"
                )
            if re.search(r"[^AGCT]", l):
                raise DREEMInputException(
                    f"reference sequences must contain only AGCT characters."
                    f" line {i} is not consisetnt with this in fasta"
                )
    log.info(f"found {num} valid reference sequences in {fa}")


def validate_csv_file(fa: str, csv: str) -> None:
    """
    ensure that the fasta file is in the correct format
    :param fa: path to the fasta file
    """
    ref_seqs = fasta_to_dict(fa)
    df = pd.read_csv(csv)
    if "name" not in df.columns:
        raise DREEMInputException(
            f"csv file does not contain a column named 'name'."
        )
    if "sequence" not in df.columns:
        raise DREEMInputException(
            f"csv file does not contain a column named 'sequence'"
        )
    if "structure" not in df.columns:
        raise DREEMInputException(
            f"csv file does not contain a column named 'structure'"
        )
    if len(ref_seqs) != len(df):
        raise DREEMInputException(
            f"number of reference sequences in fasta ({len(ref_seqs)}) does not match"
            f" number of reference sequences in dot-bracket file ({len(df)})"
        )
    for i, row in df.iterrows():
        if row["name"] not in ref_seqs:
            raise DREEMInputException(
                f"reference sequence name {row['name']} in csv file does not match"
                f" any reference sequence names in fasta file"
            )


def validate_inputs(fa, fq1, fq2, csv) -> Inputs:
    if not os.path.isfile(fa):
        raise DREEMInputException(f"fasta file: does not exist {fa}!")
    else:
        log.info(f"fasta file: {fa} exists")
        validate_fasta_file(fa)
    if not os.path.isfile(fq1):
        raise DREEMInputException(f"fastq1 file: does not exist {fq1}!")
    else:
        log.info(f"fastq1 file: {fq1} exists")
    if fq2 != "":
        if not os.path.isfile(fq2):
            raise DREEMInputException(f"fastq2 file: does not exist {fq2}!")
        else:
            log.info(f"fastq2 file: {fq2} exists")
            log.info("two fastq files supplied, thus assuming paired reads")
    if csv != "":
        if not os.path.isfile(csv):
            raise DREEMInputException(f"csv file: does not exist {csv}!")
        else:
            log.info(f"csv file: {csv} exists")
            validate_csv_file(fa, csv)
    return Inputs(fa, fq1, fq2, csv)

# run #########################################################################

def run(fasta, fastq1, fastq2, dot_bracket, params=None):
    ins = validate_inputs(fasta, fastq1, fastq2, dot_bracket)
    if params is None:
        params = get_default_params()
    else:
        validate_parameters(params)
    # perform read mapping to reference sequences
    m = Mapper()
    m.setup(params)
    m.check_program_versions()
    m.run(ins)
    # convert aligned reads to bit vectors
    bt = BitVectorGenerator()
    bt.setup(params)
    sam_path = os.path.join(
        params["dirs"]["output"], "Mapping_Files", "aligned.sam"
    )
    bt.run(sam_path, ins.fasta, ins.is_paired(), ins.csv)
    # log parameter file
    with open(os.path.join(params["dirs"]["log"], "params.yml"), "w") as f:
        yaml.dump(params, f)
