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

# TODO make sure to validate files
# TODO make sure to validate right format of fasta file "> " vs ">"
# TODO check the version of trim_galore


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

    args = parser.parse_args()
    return args


def write_log_file(fname, output):
    pass


# pipeline running functions for mapping ##############################################


def run_trim_glore(p: Parameters):
    if os.path.isfile(p.files.tg_fastq1):
        log.info(
            "SKIPPING trim_glore, it has been run already! specify -overwrite to rerun"
        )
        return

    if p.map.skip_trim_galore:
        log.info(
            "SKIPPING trim_glore, was requested by user using param map.skip_trim_glore"
        )
        log.info("copying input fastq files directly")
        shutil.copy(p.ins.fastq1, p.files.tg_fastq1)
        if p.paired:
            shutil.copy(p.ins.fastq2, p.files.tg_fastq2)
        return

    if p.paired:
        trim_galore_cmd = "trim_galore --fastqc --paired {fq1} {fq2} -o {out}".format(
            out=p.dirs.mapping, fq1=p.ins.fastq1, fq2=p.ins.fastq2
        )
    else:
        trim_galore_cmd = "trim_galore --fastqc {fq1} -o {out}".format(
            out=p.dirs.mapping, fq1=p.ins.fastq1
        )
    log.info("running trim_galore")
    log.debug(trim_galore_cmd)
    output, error = run_command(trim_galore_cmd)
    log_error("trim_galore", error)


def run_bowtie_build(p: Parameters):
    # TODO assume I already checked is a valid fasta file format
    if os.path.isfile(p.dirs.input + p.ins.ref_fasta_name + ".1.bt2"):
        log.info(
            "SKIPPING bowtie-build, it has been run already! specify -overwrite to rerun"
        )
        return

    bowtie_build_cmd = "bowtie2-build {fa} input/{fa_name}".format(
        fa=p.ins.ref_fasta, fa_name=p.ins.ref_fasta_name
    )
    log.info("running bowtie2 build")
    log.debug(bowtie_build_cmd)
    output, error = run_command(bowtie_build_cmd)
    log_error(log, "bowtie2-build", error)


# TODO valdiate bt2 args?
def run_bowtie_alignment(p: Parameters):
    if os.path.isfile(p.files.bt2_alignment_output):
        log.info(
            "SKIPPING bowtie alignment, it has been run already! specify -overwrite to rerun"
        )
        return
    bt2_args = " ".join(p.map.bt2_alignment_args.split(","))
    if p.paired:
        bowtie2_cmd = (
            "bowtie2 " + bt2_args + " -x {btindex} -1 {fq1} -2 {fq2} -S {samfile}"
        ).format(
            btindex=p.files.bt2_index,
            fq1=p.files.tg_fastq1,
            fq2=p.files.tg_fastq2,
            samfile=p.files.bt2_alignment_output,
        )
    else:
        bowtie2_cmd = (
            "bowtie2 " + bt2_args + " -x {btindex} -U {fq} -S {samfile}"
        ).format(
            btindex=p.files.bt2_index,
            fq=p.files.tg_fastq1,
            samfile=p.files.bt2_alignment_output,
        )
    log.info("running bowtie2 alignment")
    log.debug(bowtie2_cmd)
    output, error = run_command(bowtie2_cmd)
    log_error(log, "bowtie2", error)
    output_lines = output.split("\n")
    keep = []
    for l in output_lines:
        if len(l) == 0:
            continue
        if l[0] != "U":
            keep.append(l)
    log.info("results for bowtie alignment: \n" + "\n".join(keep))


def run_picard_bam_convert(p: Parameters):
    print("Converting BAM file to SAM file format")
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
    m = mapping.Mapper()
    m.run(p)
    exit()

    run_fastqc(p)
    run_trim_glore(p)
    run_bowtie_build(p)
    run_bowtie_alignment(p)
    run_picard_bam_convert(p)
    exit()

    # bit_vector.get_bit_vectors(fl, p)


if __name__ == "__main__":
    main()
