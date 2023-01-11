import sys
import click
from click_option_group import optgroup

from dreem import mapping
from dreem.logger import setup_applevel_logger, get_logger
from dreem.parameters import Inputs
from dreem.util import *


# TODO make sure to validate right format of fasta file "> " vs ">"
# TODO check the version of trim_galore
# TODO add setup.py which checks for everything

# logging/settings/args ###############################################################

log = get_logger("RUN")

def write_log_file(fname, output):
    pass


# universal logger
log = logger.get_logger("MAIN")


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
        **args
    )


def run(fasta, fastq1, fastq2, dot_bracket, parameter_file, **kwargs):
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    ins = Inputs(fasta, fastq1, fastq2, dot_bracket, parameter_file)
    #log.setLevel(logger.str_to_log_level(kwargs["log_level"]))
    # setup parameters
    # setup_parameters(args)
    # p = get_parameters()
    build_directories()
    #p.to_yaml_file(p.dirs.log + "/parameters.yml")
    # perform read mapping to reference sequences
    m = mapping.Mapper()
    m.check_program_versions()
    m.run(ins)
    #m.run(p)
    # convert aligned reads to bit vectors
    #bt = bit_vector.BitVectorGenerator()
    #bt.run(p)


if __name__ == "__main__":
    main()
