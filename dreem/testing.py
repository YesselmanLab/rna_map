"""
functions for testing dreem code
"""
import os
import shutil

from dreem.settings import get_test_path
from dreem.parameters import Inputs


def setup_directories(cur_dir):
    """
    Set up the directory for testing
    """
    os.makedirs(os.path.join(cur_dir, "input"))
    os.makedirs(os.path.join(cur_dir, "log"))
    os.makedirs(os.path.join(cur_dir, "output"))
    os.makedirs(os.path.join(cur_dir, "output", "Mapping_Files"))
    os.makedirs(os.path.join(cur_dir, "output", "BitVector_Files"))


def remove_directories(cur_dir):
    """
    Remove the directory for testing
    """
    shutil.rmtree(os.path.join(cur_dir, "input"))
    shutil.rmtree(os.path.join(cur_dir, "log"))
    shutil.rmtree(os.path.join(cur_dir, "output"))


def get_test_inputs_paired():
    """
    Get the test inputs
    """
    test_data_dir = os.path.join(get_test_path(), "resources", "case_1")
    return Inputs(
        test_data_dir + "/test.fasta",
        test_data_dir + "/test_mate1.fastq",
        test_data_dir + "/test_mate2.fastq",
    )

