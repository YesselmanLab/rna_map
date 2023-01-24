"""
test mapping functions
"""
import os
import pytest
import shutil

from rna_map.logger import setup_applevel_logger
from rna_map.mapping import Mapper
from rna_map.parameters import Inputs, get_default_params

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def get_test_inputs_paired():
    """
    Get the test inputs
    """
    test_data_dir = os.path.join(TEST_DIR, "resources", "case_1")
    return {
        "fasta": test_data_dir + "/test.fasta",
        "fastq1": test_data_dir + "/test_mate1.fastq",
        "fastq2": test_data_dir + "/test_mate2.fastq",
    }


def remove_directories(cur_dir):
    """
    Remove the directory for testing
    """
    shutil.rmtree(os.path.join(cur_dir, "input"))
    shutil.rmtree(os.path.join(cur_dir, "log"))
    shutil.rmtree(os.path.join(cur_dir, "output"))


def test_check_program_versions():
    """
    test check_program_versions
    """
    m = Mapper()
    try:
        m.check_program_versions()
    except Exception as e:
        pytest.fail(str(e))


def test_mapping():
    """
    test mapping
    """
    p = get_test_inputs_paired()
    ins = Inputs(p["fasta"], p["fastq1"], p["fastq2"], "")
    params = get_default_params()
    m = Mapper()
    m.setup(params)
    m.run(ins)
    assert os.path.isfile("input/test.1.bt2")
    try:
        m.run(ins)
    except Exception as e:
        pytest.fail(str(e))
    remove_directories(os.getcwd())
