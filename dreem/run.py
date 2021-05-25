import sys
import logging
import colorlog
import os
import subprocess
import shutil
import click

from dreem import settings, logger, mapping, bit_vector
from dreem.parameters import *
from dreem.util import *


# TODO make sure to validate right format of fasta file "> " vs ">"
# TODO check the version of trim_galore
# TODO add setup.py which checks for everything

# logging/settings/args ###############################################################


def write_log_file(fname, output):
    pass


# universal logger
log = logger.init_logger("DREEM")


def build_directories(p: Parameters):
    log.info("building directory structure")
    safe_mkdir(p.dirs.input)
    safe_mkdir(p.dirs.log)
    safe_mkdir(p.dirs.output)
    safe_mkdir(p.dirs.mapping)
    safe_mkdir(p.dirs.bitvector)


@click.command()
# required options
@click.option(
        "-fa",
        "--fasta",
        type=click.Path(exists=True),
        required=True,
        help="reference sequences in fasta format",
)
@click.option(
        "-fq1",
        "--fastq1",
        type=click.Path(exists=True),
        required=True,
        help="fastq sequencing file of mate 1",
)
@click.option(
        "-fq2",
        "--fastq2",
        type=click.Path(exists=True),
        required=False,
        help="fastq sequencing file of mate 1",
)
@click.option(
        "-ow",
        "--overwrite",
        is_flag=True,
        help="overwrites previous results, if not set will keep previous calculation "
             + "checkpoints",
)
@click.option(
        "--csv",
        required=False)
@click.option(
        "-pf", "--param-file", type=click.Path(exists=True), help="parameter files cans upply ",
)
@click.option(
        "-ll", "--log-level", help="set log level", default="INFO"
)
def main(**args):
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    log.setLevel(str_to_log_level(args["log_level"]))
    # setup parameters
    setup_parameters(args)
    p = get_parameters()
    build_directories(p)
    # perform read mapping to reference sequences
    m = mapping.Mapper()
    m.run(p)
    # convert aligned reads to bit vectors
    bt = bit_vector.BitVectorGenerator()
    bt.run(p)



if __name__ == "__main__":
    main()
