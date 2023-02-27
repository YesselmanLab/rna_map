"""
run main program of rna_map
"""

import re
import os
import pandas as pd
import yaml
from pathlib import Path

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


def validate_fasta_file(fa: str) -> bool:
    """
    ensure that the fasta file is in the correct format
    :param fa: path to the fasta file
    :return: True if valid
    """
    with open(fa, "r", encoding="utf8") as f:
        lines = f.readlines()
    if len(lines) > 2000:
        log.warning(
            "fasta file contains more than 1000 sequences, this may lead"
            " to file generation issues. Its recommended to use --summary-output-only "
        )
    num = 0
    for i, line in enumerate(lines):
        line = line.rstrip()
        if len(line) == 0:
            raise DREEMInputException(
                f"blank line found on ln: {i}. These are not allowed in fastas."
            )
        # should be a reference sequence declartion
        if i % 2 == 0:
            num += 1
            if not line.startswith(">"):
                raise DREEMInputException(
                    "reference sequence names are on line zero and even numbers."
                    f" line {i} has value which is not correct format in the fasta"
                )
            if line.startswith("> "):
                raise DREEMInputException(
                    "there should be no spaces between > and reference name."
                    f"this occured on ln: {i} in the fasta file"
                )
        elif i % 2 == 1:
            if line.startswith(">"):
                raise DREEMInputException(
                    "sequences should be on are on odd line numbers."
                    f" line {i} has value which is not correct format in fasta file"
                )
            if re.search(r"[^AGCT]", line):
                raise DREEMInputException(
                    "reference sequences must contain only AGCT characters."
                    f" line {i} is not consisetnt with this in fasta"
                )
    log.info(f"found {num} valid reference sequences in {fa}")
    return True


def validate_csv_file(fa: str, csv: str) -> None:
    """
    ensure that the fasta file is in the correct format
    :param fa: path to the fasta file
    """
    ref_seqs = fasta_to_dict(fa)
    df = pd.read_csv(csv)
    if "name" not in df.columns:
        raise DREEMInputException("csv file does not contain a column named 'name'.")
    if "sequence" not in df.columns:
        raise DREEMInputException("csv file does not contain a column named 'sequence'")
    if "structure" not in df.columns:
        raise DREEMInputException(
            "csv file does not contain a column named 'structure'"
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
                " any reference sequence names in fasta file"
            )


def validate_fastq_file(fastq_file: Path) -> bool:
    """
    validate a fastq file
    """
    #TODO add validation for gz files
    if fastq_file.suffix == ".gz":
        return True
    with open(fastq_file, "r") as f:
        lines = [f.readline().strip() for _ in range(4)]
    if len(lines) < 4:
        return False
    if not lines[0].startswith("@"):
        return False
    if not lines[2].startswith("+"):
        return False
    if len(lines[1]) != len(lines[3]):
        return False
    return True


def validate_inputs(fa, fq1, fq2, csv) -> Inputs:
    """
    validate the input files
    :param fa: path to the fasta file
    :param fq1: path to the first fastq file
    :param fq2: path to the second fastq file
    :param csv: path to the dot-bracket file
    """
    if not fa.is_file():
        raise DREEMInputException(f"fasta file: does not exist {fa}!")
    else:
        log.info(f"fasta file: {fa} exists")
        if validate_fasta_file(fa):
            log.info("fasta file is valid")
    if not fq1.is_file():
        raise DREEMInputException(f"fastq1 file: does not exist {fq1}!")
    else:
        if not validate_fastq_file(fq1):
            raise DREEMInputException(f"fastq1 file: is not a valid fastq file {fq1}!")
        log.info(f"fastq1 file: {fq1} exists")
    if str(fq2) != ".":
        if not os.path.isfile(fq2):
            raise DREEMInputException(f"fastq2 file: does not exist {fq2}!")
        else:
            log.info(f"fastq2 file: {fq2} exists")
            log.info("two fastq files supplied, thus assuming paired reads")
            if not validate_fastq_file(fq2):
                raise DREEMInputException(
                    f"fastq2 file: is not a valid fastq file {fq2}!"
                )
    if str(csv) != ".":
        if not os.path.isfile(csv):
            raise DREEMInputException(f"csv file: does not exist {csv}!")
        else:
            log.info(f"csv file: {csv} exists")
            validate_csv_file(fa, csv)
    return Inputs(fa, fq1, fq2, csv)


# run #########################################################################


def run(fasta, fastq1, fastq2, dot_bracket, params=None):
    """
    run the pipeline
    :param fasta: path to the fasta file
    :param fastq1: path to the first fastq file
    :param fastq2: path to the second fastq file
    :param dot_bracket: path to the dot-bracket file
    :param params: dictionary of parameters
    """
    ins = validate_inputs(Path(fasta), Path(fastq1), Path(fastq2), Path(dot_bracket))
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
    sam_path = Path(params["dirs"]["output"]) / "Mapping_Files" / "aligned.sam"
    bt.run(sam_path, ins.fasta, ins.is_paired(), ins.csv)
    # log parameter file
    with open(Path(params["dirs"]["log"]) / "params.yml", "w") as f:
        yaml.dump(params, f)
