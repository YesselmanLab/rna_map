"""
test running external commands
"""
import os
import pytest

from dreem.testing import (
    setup_directories,
    remove_directories,
    get_test_inputs_paired,
)
from dreem.external_cmd import (
    run_command,
    run_fastqc,
    run_trim_glore,
    run_bowtie_build,
    run_bowtie_alignment,
    run_picard_bam_convert,
    validate_bowtie2_args,
)
from dreem.exception import DREEMInputException
from dreem.logger import setup_applevel_logger

def test_run_command():
    """
    test run_command
    """
    cmd = "ls -l"
    output, error_msg = run_command(cmd)
    print(output, error_msg)

def _test_run_fastqc():
    """
    test run_fastqc
    """
    setup_directories(os.getcwd())
    ins = get_test_inputs_paired()
    run_fastqc(ins.fastq1, ins.fastq2, "output/Mapping_Files/")
    remove_directories(os.getcwd())


def _test_run_trim_glore():
    """
    test run_trim_glore
    """
    log = setup_applevel_logger()
    log.setLevel("DEBUG")
    setup_directories(os.getcwd())
    ins = get_test_inputs_paired()
    run_trim_glore(ins.fastq1, ins.fastq2, "output/Mapping_Files/")
    remove_directories(os.getcwd())


def _test_run_bowtie_build():
    """
    test run_trim_glore
    """
    setup_directories(os.getcwd())
    ins = get_test_inputs_paired()
    run_bowtie_build(ins.fasta, "input")
    assert os.path.isfile("input/test.1.bt2")
    remove_directories(os.getcwd())


def test_validate_bt2_args():
    """
    test validate_bt2_args
    """
    args = "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-p 16"
    try:
        validate_bowtie2_args(args)
    except DREEMInputException as exc:
        pytest.fail(exc)


# TODO make sure that I check for core options that should not be included
# such as -x -1 -2 -S
def test_validate_bt2_args_exceptions():
    """
    test validate_bt2_args exceptions
    """
    args = ""
    validate_bowtie2_args(args)
    # invalid boolean argument
    args = "--fake,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-p 16"
    with pytest.raises(DREEMInputException):
        validate_bowtie2_args(args)
    # invalid integer argument
    args = "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-fake 16"
    with pytest.raises(DREEMInputException):
        validate_bowtie2_args(args)
    # invalid argument type
    args = "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-p test"
    with pytest.raises(DREEMInputException):
        validate_bowtie2_args(args)


def _test_bowtie_alignment():
    """
    test bowtie alignment
    """
    log = setup_applevel_logger()
    log.setLevel("DEBUG")
    setup_directories(os.getcwd())
    ins = get_test_inputs_paired()
    run_bowtie_build(ins.fasta, "input")
    args = "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-p 16"
    run_bowtie_alignment(
        ins.fasta, ins.fastq1, ins.fastq2, "input", "output/Mapping_Files", args
    )
    assert os.path.isfile("output/Mapping_Files/aligned.sam")
    assert os.path.getsize("output/Mapping_Files/aligned.sam") > 1000
    remove_directories(os.getcwd())
