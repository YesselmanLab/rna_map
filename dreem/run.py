import re
import sys
import click
from click_option_group import optgroup

from dreem import mapping
from dreem.settings import get_py_path

from dreem.exception import DREEMInputException
from dreem.logger import setup_applevel_logger, get_logger
from dreem.parameters import Inputs, parse_parameters_from_file
from dreem.util import *

log = get_logger("RUN")


# TODO make sure to validate right format of fasta file "> " vs ">"

# logging/settings/args ###############################################################


def validate_fasta_file(fa):
    f = open(fa)
    lines = f.readlines()
    f.close()
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


def build_directories():
    log.info("building directory structure")
    os.makedirs("input")
    os.makedirs("log")
    os.makedirs("output")
    os.makedirs(os.path.join("output", "Mapping_Files"))
    os.makedirs(os.path.join("output", "BitVector_Files"))


@click.command()
@optgroup.group("main arguments")
@optgroup.option(
    "-fa",
    "--fasta",
    type=click.Path(exists=True),
    required=True,
    help="reference sequences in fasta format",
)
@optgroup.option(
    "-fq1",
    "--fastq1",
    type=click.Path(exists=True),
    required=True,
    help="fastq sequencing file of mate 1",
)
@optgroup.option(
    "-fq2",
    "--fastq2",
    type=click.Path(exists=True),
    help="fastq sequencing file of mate 2",
    default=None,
)
@optgroup.group("common options")
@optgroup.option(
    "--dot_bracket",
    type=click.Path(exists=True),
    help="A csv formatted file that contains dot bracket info for each sequence",
)
@optgroup.option(
    "-pf",
    "--param-file",
    type=click.Path(exists=True),
    help="A yml formatted file to specify parameters",
)
def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    setup_applevel_logger()
    run(
        args.pop("fasta"),
        args.pop("fastq1"),
        args.pop("fastq2"),
        args.pop("dot_bracket"),
        args.pop("param_file"),
        **args,
    )


def run(fasta, fastq1, fastq2, dot_bracket, parameter_file, **kwargs):
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    ins = Inputs(fasta, fastq1, fastq2, dot_bracket)
    # params =
    # log.setLevel(logger.str_to_log_level(kwargs["log_level"]))
    # setup parameters
    # setup_parameters(args)
    # p = get_parameters()
    build_directories()
    # p.to_yaml_file(p.dirs.log + "/parameters.yml")
    # perform read mapping to reference sequences
    m = mapping.Mapper()
    m.check_program_versions()
    m.run(ins)
    # m.run(p)
    # convert aligned reads to bit vectors
    # bt = bit_vector.BitVectorGenerator()
    # bt.run(p)


if __name__ == "__main__":
    main()
