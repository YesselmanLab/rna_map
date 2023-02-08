import os
import shutil
import pytest
from pathlib import Path

from rna_map.run import validate_inputs, validate_fasta_file
from rna_map.exception import DREEMInputException

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def get_test_inputs_paired():
    """
    Get the test inputs
    """
    test_data_dir = Path(TEST_DIR) / "resources" / "case_1"
    return {
        "fasta": test_data_dir / "test.fasta",
        "fastq1": test_data_dir / "test_mate1.fastq",
        "fastq2": test_data_dir / "test_mate2.fastq",
    }


def remove_directories(cur_dir):
    """
    Remove the directory for testing
    """
    shutil.rmtree(os.path.join(cur_dir, "input"))
    shutil.rmtree(os.path.join(cur_dir, "log"))
    shutil.rmtree(os.path.join(cur_dir, "output"))


def test_input_validation():
    """
    test input validation
    """
    p = get_test_inputs_paired()
    ins = validate_inputs(p["fasta"], p["fastq1"], "", "")
    assert p["fasta"] == ins.fasta
    assert ins.csv == ""
    assert ins.is_paired() == False
    assert ins.supplied_csv() == False

    # check to make sure we get the proper errors for supplying file that does
    # not exist
    with pytest.raises(DREEMInputException):
        validate_inputs(p["fasta"], "", "", "")
    with pytest.raises(DREEMInputException):
        validate_inputs("fake_path", p["fastq1"], "", "")
    with pytest.raises(DREEMInputException):
        validate_inputs(p["fasta"], p["fastq1"], "fake_path", "")
    with pytest.raises(DREEMInputException):
        validate_inputs(p["fasta"], p["fastq1"], "", "fake_path")


def test_fasta_checks():
    fasta_test_path = Path(TEST_DIR) / "resources" / "test_fastas"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_fasta_file(fasta_test_path / "blank_line.fasta")
    assert (
        exc_info.value.args[0]
        == "blank line found on ln: 1. These are not allowed in fastas."
    )
    with pytest.raises(DREEMInputException) as exc_info:
        validate_fasta_file(fasta_test_path / "incorrect_format.fasta")
    assert (
        exc_info.value.args[0]
        == "reference sequence names are on line zero and even numbers. line 0 "
        "has value which is not correct format in the fasta"
    )
    with pytest.raises(DREEMInputException):
        validate_fasta_file(fasta_test_path / "incorrect_sequence.fasta")
    with pytest.raises(DREEMInputException):
        validate_fasta_file(fasta_test_path / "is_rna.fasta")
    with pytest.raises(DREEMInputException):
        validate_fasta_file(fasta_test_path / "space_between_carrot.fasta")