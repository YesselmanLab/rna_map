"""
module to run commands to map sequencing reads to reference sequences
"""

import os

from pathlib import Path
import pandas as pd
from dreem import logger
from dreem.parameters import Inputs
from dreem.settings import get_py_path
from dreem.exception import (
    DREEMInputException,
    DREEMMissingRequirementsException,
    DREEMExternalProgramException,
)
from dreem.util import (
    does_program_exist,
    run_command,
    get_fastqc_version,
    get_bowtie2_version,
    get_cutadapt_version,
    get_trim_galore_version,
)

log = logger.get_logger("MAPPING")


def skip_without_overwrite(method_name):
    """
    Skip a method if the output directory already exists and overwrite is False
    """
    log.info(
        f"SKIPPING {method_name}, it has been run already! specify -overwrite "
        f"to rerun"
    )


def skip_method_by_user(method_name, method_param):
    """
    Skip a method if the user requests it
    """
    log.info(
        f"SKIPPING {method_name}, was requested by user using param "
        f"{method_param}"
    )


def run_mapping_command(method_name, cmd):
    """
    Run a mapping command and log the output
    """
    log.info(f"running {method_name}")
    log.debug(cmd)
    output, error_msg = run_command(cmd)
    if error_msg is not None:
        log.error(f"error running command: {method_name}")
        raise DREEMExternalProgramException(error_msg)
    log.info(f"{method_name} ran without errors")
    return output, error_msg


def run_fastqc(fastq1, fastq2, out_dir, overwrite=False):
    """
    run fastqc appliction on fastq files
    """
    fastqc_dir = os.path.join(out_dir, "fastqc")
    if os.path.isdir(fastqc_dir) and not overwrite:
        skip_without_overwrite("fastqc")
        return
    os.makedirs(fastqc_dir)
    fastqc_cmd = f"fastqc {fastq1} {fastq2} -o {fastqc_dir}"
    run_mapping_command("fastqc", fastqc_cmd)


def run_trim_glore(fastq1, fastq2, out_dir):
    """
    Run trim glore on fastq files
    """
    if fastq2 == "":
        cmd = f"trim_galore --fastqc --paired {fastq1} {fastq2} -o {out_dir}"
    else:
        cmd = f"trim_galore --fastqc {fastq1} -o {out_dir}"
    run_mapping_command("trim_galore", cmd)


def run_bowtie_build(fasta, input_dir, overwrite=False):
    fasta_name = Path(fasta).stem
    bt2_file = os.path.join(input_dir, fasta_name + ".1.bt2")
    if os.path.isfile(bt2_file) and not overwrite:
        skip_without_overwrite("bowtie-build")
        return
    cmd = f'bowtie2-build "{fasta}" input/{fasta_name}'
    run_mapping_command("bowtie2-build", cmd)


def validate_bowtie2_args(args: str):
    """
    Validate the bowtie2 arguments
    """

    def check_type(arg):
        """
        Check the type of an argument
        :arg: the argument to check
        :return: the type of the argument
        """
        if arg.isdigit():
            return "int"
        try:
            float(arg)
            return "float"
        except ValueError:
            return "str"

    df = pd.read_csv(get_py_path() + "resources/bowtie2_args.csv")
    valid_bt2_args = {}
    for _, row in df.iterrows():
        valid_bt2_args[row["param"]] = row["vtype"]
    if len(args) == 0:
        log.warning("no bowtie2 arguments supplied thats probably wrong")
    supplied_args = args.strip().split(",")
    for full_arg in supplied_args:
        if len(full_arg) == 0:
            continue
        if full_arg in valid_bt2_args:
            log.debug(f"{full_arg} is a valid bt2 argument")
            continue
        spl = full_arg.split()
        if len(spl) == 1:
            raise DREEMInputException(
                f"{full_arg} is not a valid bowtie2 argument. "
                f"Please check the documentation for valid arguments"
            )
        arg, arg_val = spl[0], spl[1]
        if arg in valid_bt2_args:
            log.debug(f"{arg} is a valid bt2 argument")
        else:
            raise DREEMInputException(f"{full_arg} is an invalid bt2 argument")
        if check_type(arg_val) != valid_bt2_args[arg]:
            raise DREEMInputException(
                f"{arg} must be of type {valid_bt2_args[arg]}"
            )
    log.debug("all bt2 arguments are valid")


def run_bowtie_alignment(fasta, fastq1, fastq2, in_dir, out_dir, args):
    """
    Run bowtie2 alignment
    """
    # check to make sure bt2 args are valid
    validate_bowtie2_args(args)
    bt2_index = in_dir + "/" + Path(fasta).stem
    bt2_args = " ".join(args.split(","))
    sam_file = out_dir + "/aligned.sam"
    cmd = f"bowtie2 {bt2_args} -x {bt2_index} -S {sam_file} "
    if fastq2 != "":
        cmd += f"-1 {fastq1} -2 {fastq2}"
    else:
        cmd += f"-U {fastq1}"
    output, _ = run_mapping_command("bowtie2 alignment", cmd)
    output_lines = output.split("\n")
    keep = []
    for l in output_lines:
        if len(l) == 0:
            continue
        if l[0] != "U":
            keep.append(l)
    log.info("results for bowtie alignment: \n" + "\n".join(keep))


def run_picard_bam_convert(sam_file, bam_file):
    log.info("Converting BAM file to SAM file format")
    picard_path = get_py_path() + "/resources/picard.jar"
    cmd = (
        f"java -jar {picard_path} SamFormatConverter I={sam_file} O={bam_file}"
    )
    run_mapping_command("picard BAM conversion", cmd)


def run_picard_sort(bam_file, sorted_bam_file):
    log.info("sorting BAM file")
    picard_path = get_py_path() + "/resources/picard.jar"
    cmd = (
        f"java -jar {picard_path} SORT_ORDER=coordinate I={bam_file} "
        f"O={sorted_bam_file}"
    )
    run_mapping_command("picard BAM sort", cmd)


class Mapper(object):
    def check_program_versions(self):
        """
        check to make sure all required programs exist and have the
        right versions
        """
        programs = "fastqc,bowtie2,trim_galore,cutadapt".split(",")
        for prog in programs:
            if not does_program_exist(prog):
                self.__program_not_found(prog)

        log.info(f"bowtie2 {get_bowtie2_version()} detected!")
        log.info(f"fastqc {get_fastqc_version()} detected!")
        log.info(f"trim_galore {get_trim_galore_version()} detected!")
        log.info(f"cutapt {get_cutadapt_version()} detected!")

    def run(self, ins: Inputs, params):
        pass
        # don't rerun unless asked with -overwrite
        # if os.path.isfile(self._p.files.tg_fastq1) and not self._p.map.overwrite:
        #    self.__skip_without_overwrite("trim_galore")
        #    return
        #run_fastqc(ins.fastq1, ins.fastq2, out_dir + "/Mapping_Files/")

        # skip by user request, make sure to copy the files for next step
        # skip trim galore
        """
        if skip:
            skip_method_by_user("trim_glore", "skip_trim_galore")
            log.info("copying input fastq files directly")
            shutil.copy(self._p.ins.fastq1, self._p.files.tg_fastq1)
            if self._p.paired:
                shutil.copy(self._p.ins.fastq2, self._p.files.tg_fastq2)
            return
        """
        # self.__run_trim_glore()  # run trim galore
        # self.__run_bowtie_build()  # run bowtie build
        # self.__run_bowtie_alignment()  # row bowtie
        # self.__run_picard_bam_convert()  # convert sam to bam
        # self.__run_picard_sort()
        # self.__run_picard_metrics()
        log.info("finished mapping!")

    def __program_not_found(self, p_name):
        msg = f"{p_name} is not found make sure it is accessible in $PATH"
        raise DREEMMissingRequirementsException(msg)

    def __skip_method_by_user(self, method_name, method_param):
        log.info(
            "SKIPPING {}, was requested by user using param {}".format(
                method_name, method_param
            )
        )

    def __log_output(self):
        # f = open(f"{self._p.dirs.log}/{method_name}.log", "w")
        # f.write(output)
        # f.close()
        pass

    # run programs #############################################################
