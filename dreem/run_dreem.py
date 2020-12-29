import sys
import argparse
import logging
import colorlog
import os
import subprocess

from dreem import BitVector, settings, logger
from dreem.parameters import *

# TODO make sure to validate files
# TODO make sure to validate right format of fasta file "> " vs ">"
# TODO check the version of trim_galore


# logging/settings/args ###############################################################


def safe_mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


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

    args = parser.parse_args()
    return args


def run_command(cmd):
    output, error_msg = None, None
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT
        ).decode("utf8")
    except subprocess.CalledProcessError as exc:
        error_msg = exc.output.decode("utf8")
    return output, error_msg


def write_log_file(fname, output):
    pass


def log_error(pname, error_msg):
    if error_msg is None:
        log.info("{} ran without errors".format(pname))
    else:
        log.error("{}} returned error:".format(pname))
        log.error(error_msg)
        log.error("EXITING")
        exit()


# pipeline running functions for mapping ##############################################


def run_fastqc(p: Parameters):
    # run fastqc
    fastqc_dir = p.dirs.mapping + "/" + p.ins.fastq1_name + "_fastqc"
    if os.path.isdir(fastqc_dir) and not p.overwrite:
        log.info(
            "SKIPPING fastqc, it has been run already! specify -overwrite to rerun"
        )
        return
    if p.map.skip_fastqc:
        log.info("SKIPPING fastqc, was requested by user using param map.skip_fastqc")
        return

    if p.paired:
        fastqc_cmd = "fastqc --extract {fq1} {fq2} --outdir={dir}".format(
            dir=p.dirs.mapping, fq1=p.ins.fastq1, fq2=p.ins.fastq2,
        )
    else:
        fastqc_cmd = "fastqc --extract {fq1} -outdir={dir}".format(
            dir=p.dirs.mapping, fq1=p.ins.fastq1
        )

    log.info("running fastqc")
    log.debug(fastqc_cmd)
    output, error = run_command(fastqc_cmd)
    log_error("fastqc", error)


def run_trim_glore(p: Parameters):
    if os.path.isfile(p.files.tg_fastq1):
        log.info(
            "SKIPPING trim_glore, it has been run already! specify -overwrite to rerun"
        )
        return

    if p.paired:
        trim_galore_cmd = "trim_galore --fastqc --paired {fq1} {fq2} -o {out}".format(
            out=p.dirs.mapping, fq1=p.ins.fastq1, fq2=p.ins.fastq2
        )
    else:
        trim_galore_cmd = "trim_galore --fastqc {fq1} -o {out}".format(
            out=p.dirs.mapping, fq1=p.ins.fastq1
        )
    log.info("#######################################################################")
    log.info("running trim_galore")
    log.info("#######################################################################")
    log.debug(trim_galore_cmd)


def run_bowtie_build(args):
    # TODO assume I already checked is a valid fasta file format
    fa_name = args.fasta.split(".")[0]
    bowtie_build_cmd = "bowtie2-build {fa} input/{fa_name}".format(
        fa=args.fasta, fa_name=fa_name
    )
    log.info("running bowtie2 build")
    log.debug(bowtie_build_cmd)
    os.system(bowtie_build_cmd)


def run_bowtie_alignment(fl):
    bowtie2_cmd = (
        "bowtie2 --local --no-unal --no-discordant --no-mixed -X 1000 -L 12"
        + " -p 16 -x {btindex} -1 {fq1} -2 {fq2} -S {samfile}"
    ).format(
        btindex=fl.btindex,
        fq1=fl.trim_galore_fastq1,
        fq2=fl.trim_galore_fastq2,
        samfile=fl.sam_file,
    )
    log.info("running bowtie2 alignment")
    os.system(bowtie2_cmd)


def run_picard(fl):
    picard_path = settings.get_lib_path() + "/resources/picard.jar"
    picard_cmd = (
        "java -jar "
        + picard_path
        + " SamFormatConverter I={sam_file} O={bam_file}".format(
            sam_file=fl.sam_file, bam_file=fl.bam_file
        )
    )
    log.info("running picard")
    os.system(picard_cmd)


# universal logger
log = logger.init_logger("DREEM")


def build_directories(p: Parameters):
    log.info("building directory structure")
    safe_mkdir(p.dirs.input)
    safe_mkdir(p.dirs.log)
    safe_mkdir(p.dirs.output)
    safe_mkdir(p.dirs.mapping)
    safe_mkdir(p.dirs.bitvector)


def check_input_files(p: Parameters):
    pass


def main():
    log.setLevel(logging.DEBUG)
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    args = parse_args()

    setup_parameters(args)
    p = get_parameters()

    build_directories(p)
    run_fastqc(p)
    # run_trim_glore(p)
    exit()

    """
    run_bowtie_build(args)
    run_bowtie_alignment(fl)
    run_picard(fl)"""
    """picard_path = get_lib_path() + "/resources/picard.jar"
    print('Converting BAM file to SAM file format')
    convert_cmd = 'java -jar {} SamFormatConverter I={} O={} >> log.txt'
    convert_cmd = convert_cmd.format\
        (picard_path, fl.bam_file, fl.new_sam_file)
    os.system(convert_cmd)"""

    # bit_vector.get_bit_vectors(fl, p)


if __name__ == "__main__":
    main()
