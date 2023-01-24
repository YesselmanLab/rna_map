"""
test cli interface
"""
import os
import shutil

from click.testing import CliRunner
from rna_map.cli import cli

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

def remove_directories(cur_dir):
    """
    Remove the directory for testing
    """
    shutil.rmtree(os.path.join(cur_dir, "input"))
    shutil.rmtree(os.path.join(cur_dir, "log"))
    shutil.rmtree(os.path.join(cur_dir, "output"))



def test_cli_single():
    """
    test running the program
    """
    path = TEST_DIR + "/resources/case_unit/"
    runner = CliRunner()
    result = runner.invoke(
        cli, ["-fa", f"{path}/test.fasta", "-fq1", f"{path}/test_mate1.fastq"]
    )
    assert result.exit_code == 0
    #assert os.path.isfile(f"output/test_mate1.fastq")
    remove_directories(os.getcwd())
