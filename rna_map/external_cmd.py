"""
run and check external commannds
"""
import os
import shutil
import subprocess
from typing import Optional
from pathlib import Path
from dataclasses import dataclass
import pandas as pd

from rna_map.settings import get_py_path
from rna_map.logger import get_logger
from rna_map.exception import DREEMInputException, DREEMExternalProgramException

log = get_logger("EXTERNAL_CMD")


@dataclass(frozen=True, order=True)
class ProgOutput:
    """
    Class to store the output of an external program
    """

    output: Optional[str]
    error: Optional[str]


def does_program_exist(prog_name: str) -> bool:
    """
    Check if a program exists
    :prog_name: name of the program
    """
    if shutil.which(prog_name) is None:
        return False
    else:
        return True


def get_bowtie2_version() -> str:
    """
    Get the version of bowtie2
    :return: version of bowtie2
    """
    if not does_program_exist("bowtie2"):
        raise DREEMExternalProgramException(
                "cannot get bowtie2 version, cannot find the exe"
        )
    output = subprocess.check_output("bowtie2 --version", shell=True).decode(
            "utf8"
    )
    lines = output.split("\n")
    l_spl = lines[0].split()
    return l_spl[-1]


def get_fastqc_version() -> str:
    """
    Get the version of fastqc
    :return: version of fastqc
    """
    if not does_program_exist("fastqc"):
        raise DREEMExternalProgramException(
                "cannot get fastqc version, cannot find the exe"
        )
    out = run_command("fastqc --version")
    lines = out.output.split("\n")
    if len(lines) < 1:
        raise ValueError(
                "cannot get fastqc version, output is not valid: {}".format(
                        out.output
                )
        )
    l_spl = lines[0].split()
    return l_spl[1]


def get_trim_galore_version():
    if not does_program_exist("trim_galore"):
        raise DREEMExternalProgramException(
                "cannot get trim_galore version, cannot find the exe"
        )
    output = subprocess.check_output(
            "trim_galore --version", shell=True
    ).decode("utf8")
    lines = output.split("\n")
    if len(lines) < 4:
        raise ValueError(
                "cannot get fastqc version, output is not valid: {}".format(output)
        )
    for l in lines:
        if l.find("version") != -1:
            l_spl = l.split()
            return l_spl[-1]
    return ""


def get_cutadapt_version():
    if not does_program_exist("cutadapt"):
        raise DREEMExternalProgramException(
                "cannot get cutadapt version, cannot find the exe"
        )
    output = subprocess.check_output("cutadapt --version", shell=True).decode(
            "utf8"
    )
    return output.rstrip().lstrip()


def run_command(cmd: str) -> ProgOutput:
    """
    Run a command and return the output
    :cmd: command to run
    """
    output, error_msg = None, None
    try:
        output = subprocess.check_output(
                cmd, shell=True, stderr=subprocess.STDOUT
        ).decode("utf8")
    except subprocess.CalledProcessError as exc:
        error_msg = exc.output.decode("utf8")
    return ProgOutput(output, error_msg)


def run_named_command(method_name: str, cmd: str) -> ProgOutput:
    """
    Run a mapping command and log the output
    :method_name: name of the method
    :cmd: command to run
    :return: output of the command
    """
    log.info(f"running {method_name}")
    log.debug(cmd)
    out = run_command(cmd)
    if out.error is not None:
        log.error(f"error running command: {method_name}")
        raise DREEMExternalProgramException(out.error)
    log.info(f"{method_name} ran without errors")
    return out


def run_fastqc(fastq1: str, fastq2: str, out_dir: str) -> ProgOutput:
    """
    run fastqc appliction on fastq files
    :fastq1: path to fastq1 file
    :fastq2: path to fastq2 file
    :out_dir: path to output directory
    """
    fastqc_dir = os.path.join(out_dir, "fastqc")
    os.makedirs(fastqc_dir, exist_ok=True)
    fastqc_cmd = f"fastqc {fastq1} {fastq2} -o {fastqc_dir}"
    return run_named_command("fastqc", fastqc_cmd)


def run_trim_glore(fastq1: str, fastq2: str, out_dir: str) -> ProgOutput:
    """
    Run trim glore on fastq files
    :fastq1: path to fastq1 file
    :fastq2: path to fastq2 file
    :out_dir: path to output directory
    """
    if fastq2 != "":
        cmd = f"trim_galore --fastqc --paired {fastq1} {fastq2} -o {out_dir}"
    else:
        cmd = f"trim_galore --fastqc {fastq1} -o {out_dir}"
    return run_named_command("trim_galore", cmd)


def run_bowtie_build(fasta: str, input_dir: str) -> ProgOutput:
    """
    Run bowtie2-build on a fasta file
    :fasta: path to fasta file
    :input_dir: path to input directory
    """
    fasta_name = Path(fasta).stem
    cmd = f'bowtie2-build "{fasta}" {input_dir}/{fasta_name}'
    return run_named_command("bowtie2-build", cmd)


def validate_bowtie2_args(args: str) -> bool:
    """
    Validate the bowtie2 arguments
    :args: arguments to validate, seperated by ","
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


def run_bowtie_alignment(
        fasta: str, fastq1: str, fastq2: str, in_dir: str, out_dir: str, args: str,
        **kwargs) -> ProgOutput:
    """
    Run bowtie2 alignment
    :fasta: path to fasta file
    :fastq1: path to fastq1 file
    :fastq2: path to fastq2 file
    :in_dir: path to bowtie2 index directory
    """
    # check to make sure bt2 args are valid
    validate_bowtie2_args(args)
    bt2_index = in_dir + "/" + Path(fasta).stem
    bt2_args = " ".join(args.split(","))
    sam_file = out_dir + "/aligned.sam"
    cmd = f"bowtie2 {bt2_args} -x {bt2_index} -S {sam_file} "
    if fastq2 != "":
        cmd += f"-1 {fastq1} -2 {fastq2} "
    else:
        cmd += f"-U {fastq1}"
    if "save_unaligned" in kwargs:
        cmd += " --un-conc unaligned.fastq"
    out = run_named_command("bowtie2 alignment", cmd)
    output_lines = out.output.split("\n")
    keep = []
    for l in output_lines:
        if len(l) == 0:
            continue
        if l[0] != "U":
            keep.append(l)
    log.info("results for bowtie alignment: \n" + "\n".join(keep))
    return out
