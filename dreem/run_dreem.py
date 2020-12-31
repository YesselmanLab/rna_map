import sys
import argparse
import logging
import colorlog
import os
import subprocess
import shutil

from dreem import BitVector, settings, logger, mapping
from dreem.parameters import *
from dreem.util import *

# TODO make sure to validate right format of fasta file "> " vs ">"
# TODO check the version of trim_galore
# TODO add setup.py which checks for everything

# logging/settings/args ###############################################################


def parse_args():
    parser = argparse.ArgumentParser(description="")
    required_named = parser.add_argument_group("required args")
    required_named.add_argument(
        "-fa", "--fasta", help="reference sequences ", required=True
    )
    required_named.add_argument(
        "-fq1", "--fastq1", help="fastq sequencing file of mate 1 ", required=True
    )
    required_named.add_argument(
        "-fq2", "--fastq2", help="fastq sequencing file of mate 2 ", required=False
    )
    required_named.add_argument(
        "-o", "--overwrite", help=" ", action="store_true", required=False
    )
    required_named.add_argument("-p", "--params", help=" ", required=False)
    required_named.add_argument(
        "-ll", "--log_level", help=" ", default="info", required=False
    )

    args = parser.parse_args()
    return args


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


def main():
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    args = parse_args()
    log.setLevel(str_to_log_level(args.log_level))

    setup_parameters(args)
    p = get_parameters()
    build_directories(p)
    m = mapping.Mapper()
    m.run(p)
    exit()

    # bit_vector.get_bit_vectors(fl, p)


if __name__ == "__main__":
    main()
