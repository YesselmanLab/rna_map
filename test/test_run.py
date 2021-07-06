import pytest
import os
from click.testing import CliRunner

from dreem import run, settings, util

TEST_DIR = os.path.dirname( os.path.realpath( __file__ ) )
BASE_DIR = os.path.dirname(TEST_DIR)

@pytest.mark.integration
def _test_case_1():
    util.safe_rmdir("output")
    util.safe_rmdir("input")

    path = BASE_DIR + "/test/resources/case_1/"
    runner = CliRunner()

    args = [
        "--fasta", path + "test.fasta",
        "--fastq1", path + "test_mate1.fastq",
        "--fastq2", path + "test_mate2.fastq",
        "--overwrite",
    ]
    result = runner.invoke(run, args, prog_name='dreem-test')
    
    assert os.path.isfile("output/BitVector_Files/mttr-6-alt-h3_1_134_pop_avg.html")
    util.safe_rmdir("output")
    util.safe_rmdir("input")

