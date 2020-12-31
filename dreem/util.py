import os
import subprocess
import shutil


def safe_mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


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
    output = subprocess.check_output("bowtie2 --version", shell=True).decode("utf8")
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
    l_spl = lines[1].split()
    return l_spl[1]


def get_trim_galore_version():
    if not does_program_exist("trim_galore"):
        raise ValueError("cannot get trim_galore version, cannot find the exe")
    output = subprocess.check_output("trim_galore --version", shell=True).decode("utf8")
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
    output = subprocess.check_output("cutadapt --version", shell=True).decode("utf8")
    return output.rstrip().lstrip()
