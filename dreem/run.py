import sys
import logging
import colorlog
import os
import subprocess
import shutil
import click
from click_option_group import optgroup

from dreem import settings, logger, mapping, bit_vector, args
from dreem.parameters import *
from dreem.util import *


# TODO make sure to validate right format of fasta file "> " vs ">"
# TODO check the version of trim_galore
# TODO add setup.py which checks for everything

# logging/settings/args ###############################################################


def write_log_file(fname, output):
    pass


# universal logger
log = logger.init_logger("DREEM", log_outfile="dreem.log")


def build_directories(p: Parameters):
    log.info("building directory structure")
    safe_mkdir(p.dirs.input)
    safe_mkdir(p.dirs.log)
    safe_mkdir(p.dirs.output)
    safe_mkdir(p.dirs.mapping)
    safe_mkdir(p.dirs.bitvector)

@click.command()
@optgroup.group("main arguments and options")
def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    log.setLevel(str_to_log_level(args["log_level"]))
    exit()
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
