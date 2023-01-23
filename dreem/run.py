import re
import sys
import logging
import subprocess
import cloup
from cloup import option_group, option
from dreem import mapping
from dreem.settings import get_py_path

from dreem.cli_opts import (
    bit_vector_options,
    main_options,
    mapping_options,
    misc_options,
    parse_cli_args,
)
from dreem.exception import DREEMInputException
from dreem.logger import setup_applevel_logger, get_logger
from dreem.parameters import (
    Inputs,
    parse_parameters_from_file,
    get_default_params,
    validate_parameters,
)
from dreem.bit_vector import BitVectorGenerator
from dreem.external_cmd import does_program_exist
from dreem.util import *

log = get_logger("RUN")


# logging/settings/args ###############################################################


def validate_fasta_file(fa: str) -> None:
    """
    ensure that the fasta file is in the correct format
    :param fa: path to the fasta file
    """
    with open(fa, "r", encoding="utf8") as f:
        lines = f.readlines()
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


def validate_inputs(fa, fq1, fq2, csv):
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
            log.info("fastq2 file: {} exists".format(fq2))
            log.info("two fastq files supplied, thus assuming paired reads")
    if csv != "":
        if not os.path.isfile(fq2):
            raise DREEMInputException(
                "csv file: does not exist {}!".format(fq2)
            )
        else:
            log.info("csv file: {} exists".format(fq2))
    return Inputs(fa, fq1, fq2, csv)


def run_in_docker(args):
    if not does_program_exist("docker"):
        raise ValueError("docker is not installed")

    file_map = {
        "dot_bracket": "test.csv",
        "param_file": "test.yml",
        "fasta": "test.fasta",
        "fastq1": "test_mate1.fastq",
        "fastq2": "test_mate2.fastq",
    }
    file_args = {
        "dot_bracket": "--dot-bracket",
        "param_file": "--param-file",
        "fasta": "--fasta",
        "fastq1": "--fastq1",
        "fastq2": "--fastq2",
    }
    files = "fasta,fastq1,fastq2,dot_bracket,param_file".split(",")
    docker_cmd = f"docker run -v $(pwd):/data "
    dreem_cmd = "dreem "
    dirs = {os.getcwd(): "/data"}
    dcount = 2
    # TODO add other args into dreem_cmd
    for f in files:
        f_path = args[f]
        if f_path is None or f_path == "":
            continue
        dir_name = os.path.abspath(os.path.dirname(f_path))
        if dir_name == os.path.abspath(os.getcwd()):
            dreem_cmd += f"{file_args[f]} {file_map[f]} "
            continue
        if dir_name not in dirs:
            dirs[dir_name] = f"/data{dcount}"
            docker_cmd += f"-v {dir_name}:/{dirs[dir_name]} "
            dcount += 1
        dreem_cmd += f"{file_args[f]} {dirs[dir_name]}/{file_map[f]} "
    docker_cmd += " dreem_dev  "
    docker_cmd += dreem_cmd
    print(docker_cmd)


# cli #########################################################################

# TODO validate that the csv file is in the correct format
# TODO add options from command line into params
@cloup.command()
@main_options()
@mapping_options()
@bit_vector_options()
@misc_options()
def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    if args["debug"]:
        setup_applevel_logger(is_debug=True)
        log.info("Debug logging is on")
    else:
        setup_applevel_logger()
    # print(sys.argv)
    if args["docker"]:
        print("running in docker")
        run_in_docker(args)
        return
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    if args["param_file"] is None:
        params = get_default_params()
    else:
        params = parse_parameters_from_file(args["param_file"])
    parse_cli_args(params, args)
    # TODO add cmd line arguments into params
    run(
        args.pop("fasta"),
        args.pop("fastq1"),
        args.pop("fastq2"),
        args.pop("dot_bracket"),
        params,
    )


# TODO flag a warning about fasta files with more than 1000 sequences
# TODO add validation back in
# TODO validate that the csv file is in the correct format
def run(fasta, fastq1, fastq2, dot_bracket, params=None):
    ins = Inputs(fasta, fastq1, fastq2, dot_bracket)
    if params is None:
        params = get_default_params()
    else:
        validate_parameters(params)
    # p.to_yaml_file(p.dirs.log + "/parameters.yml")
    # perform read mapping to reference sequences
    m = mapping.Mapper()
    m.setup(params)
    m.check_program_versions()
    m.run(ins)
    # convert aligned reads to bit vectors
    bt = BitVectorGenerator()
    bt.setup(params)
    sam_path = os.path.join(
        params["dirs"]["output"], "Mapping_Files", "converted.sam"
    )
    bt.run(sam_path, ins.fasta, ins.is_paired(), ins.csv)


if __name__ == "__main__":
    main()
