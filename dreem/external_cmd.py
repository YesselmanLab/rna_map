"""
run and check external commannds
"""
import os
import shutil
import subprocess
from pathlib import Path
import pandas as pd

from dreem.settings import get_py_path
from dreem.logger import get_logger
from dreem.exception import DREEMInputException

log = get_logger("EXTERNAL_CMD")

def run_command(cmd):
    output, error_msg = None, None
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT
        ).decode("utf8")
    except subprocess.CalledProcessError as exc:
        error_msg = exc.output.decode("utf8")
    return output, error_msg


def does_program_exist(prog_name):
    if shutil.which(prog_name) is None:
        return False
    else:
        return True


def get_bowtie2_version():
    if not does_program_exist("bowtie2"):
        raise ValueError("cannot get bowtie2 version, cannot find the exe")
    output = subprocess.check_output("bowtie2 --version", shell=True).decode(
        "utf8"
    )
    lines = output.split("\n")
    l_spl = lines[0].split()
    return l_spl[-1]


def get_fastqc_version():
    if not does_program_exist("fastqc"):
        raise ValueError("cannot get fastqc version, cannot find the exe")
    (output, _) = run_command("fastqc --version")
    lines = output.split("\n")
    if len(lines) < 1:
        raise ValueError(
            "cannot get fastqc version, output is not valid: {}".format(output)
        )
    l_spl = lines[0].split()
    return l_spl[1]


def get_trim_galore_version():
    if not does_program_exist("trim_galore"):
        raise ValueError("cannot get trim_galore version, cannot find the exe")
    output = subprocess.check_output(
        "trim_galore --version", shell=True
    ).decode("utf8")
    lines = output.split("\n")
    if len(lines) < 4:
        raise ValueError(
            "cannot get fastqc version, output is not valid: {}".format(output)
        )
    l_spl = lines[3].split()
    return l_spl[1]


def get_cutadapt_version():
    if not does_program_exist("cutadapt"):
        raise ValueError("cannot get cutadapt version, cannot find the exe")
    output = subprocess.check_output("cutadapt --version", shell=True).decode(
        "utf8"
    )
    return output.rstrip().lstrip()



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


def run_fastqc(fastq1, fastq2, out_dir):
    """
    run fastqc appliction on fastq files
    """
    fastqc_dir = os.path.join(out_dir, "fastqc")
    os.makedirs(fastqc_dir, exist_ok=True)
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


def run_bowtie_build(fasta, input_dir):
    fasta_name = Path(fasta).stem
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