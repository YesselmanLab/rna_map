import pytest
import os
from click.testing import CliRunner

from dreem import run, settings, util, run_docker

@pytest.mark.integration
def test_case_1():
    path = settings.get_lib_path() + "/test/resources/case_1/"
    runner = CliRunner()
    args = [
        "--fasta", path + "test.fasta",
        "--fastq1", path + "test_mate1.fastq",
        "--fastq2", path + "test_mate2.fastq",
        "-ow"
    ]
    result = runner.invoke(run_docker, args, prog_name='dreem-test')
          
    assert os.path.isfile("output/BitVector_Files/mttr-6-alt-h3_1_134_pop_avg.html")

    util.safe_rmdir(path + "output")
    util.safe_rmdir(path + "input")


