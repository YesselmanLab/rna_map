import sys
import colorlog
import os
import subprocess
import shutil
import click
from click_option_group import optgroup

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
log = logger.log


def build_directories(p: Parameters):
    log.info("building directory structure")
    safe_mkdir(p.dirs.input)
    safe_mkdir(p.dirs.log)
    safe_mkdir(p.dirs.output)
    safe_mkdir(p.dirs.mapping)
    safe_mkdir(p.dirs.bitvector)


@click.command()
@optgroup.group("main arguments")
@optgroup.option("-fa", "--fasta", type=click.Path(exists=True), required=True,
                 help="reference sequences in fasta format")
@optgroup.option("-fq1", "--fastq1", type=click.Path(exists=True), required=True,
                 help="fastq sequencing file of mate 1")
@optgroup.option("-fq2", "--fastq2", type=click.Path(exists=True),
                 help="fastq sequencing file of mate 2", default=None)
@optgroup.group("common options")
@optgroup.option("--dot_bracket", type=click.Path(exists=True),
                 help="A csv formatted file that contains dot bracket info for each sequence")
@optgroup.option("-pf", "--param-file", type=click.Path(exists=True),
                 help="A yml formatted file to specify parameters")
@optgroup.option("-ow", "--overwrite", is_flag=True,
                 help="overwrites previous results, if not set will keep previous "
                      "calculation checkpoints")
@optgroup.option("-ll", "--log-level", help="set log level (INFO|WARN|DEBUG|ERROR|FATAL)",
                 default="INFO")
@optgroup.option("-rob", "--restore_org_behavior", is_flag=True, default=False,
                 help="retores the original behavior of dreem upon first release")
@optgroup.group("map options")
@optgroup.option("--map-overwrite", is_flag=True,
                 help="overwrite mapping calculation")
@optgroup.option("--skip", is_flag=True,
                 help="do not perform sequence mapping, not recommended")
@optgroup.option("--skip_fastqc", is_flag=True,
                 help="do not run fastqc for quality control of sequence data")
@optgroup.option("--skip_trim_galore", is_flag=True,
                 help="do not run trim galore to remove adapter sequences at ends")
@optgroup.option("--tg_q_cutoff", default=None,
                 help="TODO")
@optgroup.option("--bt2_alignment_args", default=None,
                 help="TODO")
@optgroup.group("bv options")
@optgroup.option("--skip", is_flag=True,
                 help="skip bit vector generation step, not recommended")
@optgroup.option("--bv-overwrite", is_flag=True,
                 help="overwrite bit vector calculation")
@optgroup.option("--qscore_cutoff", default=None,
                 help="quality score of read nucleotide, sets to ambigious if under this val")
@optgroup.option("--num_of_surbases", default=None,
                 help="number of bases around a mutation")
@optgroup.option("--map_score_cutoff", default=None,
                 help="map alignment score cutoff for a read, read is discarded if under this value")
@optgroup.option("--mutation_count_cutoff", default=None,
                 help="maximum number of mutations in a read allowable")
@optgroup.option("--percent_length_cutoff", default=None,
                 help="read is discarded if less than this percent of a ref sequence is included")
@optgroup.option("--summary_output_only", is_flag=True,
                 help="")
@optgroup.option("--plot_sequence", is_flag=True,
                 help="")
def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    run(args)


def run(args):
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    log.setLevel(logger.str_to_log_level(args["log_level"]))
    # setup parameters
    setup_parameters(args)
    p = get_parameters()
    build_directories(p)
    p.to_yaml_file(p.dirs.log + "/parameters.yml")
    # perform read mapping to reference sequences
    m = mapping.Mapper()
    m.run(p)
    # convert aligned reads to bit vectors
    bt = bit_vector.BitVectorGenerator()
    bt.run(p)


if __name__ == "__main__":
    main()
