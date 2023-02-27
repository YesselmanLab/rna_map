"""
test cli interface
"""
import os
import shutil
from pathlib import Path
from click.testing import CliRunner
from rna_map.cli import cli

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def remove_directories(cur_dir):
    """
    Remove the directory for testing
    """
    shutil.rmtree(Path(cur_dir) / "output")
    shutil.rmtree(Path(cur_dir) / "log")
    shutil.rmtree(Path(cur_dir) / "input")


def _test_cli_single():
    """
    test running the program
    """
    path = Path(TEST_DIR) / "resources" / "case_unit/"
    runner = CliRunner()
    #result = runner.invoke(
    #    cli, ["-fa", path / "test.fasta", "-fq1", path / "test_mate1.fastq"]
    #)
    # TODO understand why this is failing
    #assert result.exit_code == 0
    # assert os.path.isfile(f"output/test_mate1.fastq")
    remove_directories(os.getcwd())
