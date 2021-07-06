import pytest
import shutil
import os
from click.testing import CliRunner

from dreem import run, settings, util, run_docker

TEST_DIR = os.path.dirname( os.path.realpath( __file__ ) )
BASE_DIR = os.path.dirname(TEST_DIR)


@pytest.mark.integration
def test_case_1():
    util.safe_rmdir("output")
    util.safe_rmdir("input")

    path = BASE_DIR + "/test/resources/case_1/"
    runner = CliRunner()

    args = [
        "--fasta", path + "test.fasta",
        "--fastq1", path + "test_mate1.fastq",
        "--fastq2", path + "test_mate2.fastq",
        "-ow"
    ]
    result = runner.invoke(run_docker, args, prog_name='dreem-test')
          
    assert os.path.isfile("output/BitVector_Files/mttr-6-alt-h3_1_134_pop_avg.html")
    util.safe_rmdir("output")
    util.safe_rmdir("input")


def test_entrypoints_work():

    assert not os.system('dreem --help')
    assert not os.system('dreem-docker --help')


@pytest.mark.integration
def test_case_2():
    DIR_WITH_SPACES = 'test\\ dir'
    
    util.safe_rmdir( DIR_WITH_SPACES )
    os.mkdir(DIR_WITH_SPACES)
    os.chdir(DIR_WITH_SPACES)
    util.safe_rmdir("output")
    util.safe_rmdir("input")

    path = BASE_DIR + "/test/resources/case_1/"
    runner = CliRunner()

    args = [
        "--fasta", path + "test.fasta",
        "--fastq1", path + "test_mate1.fastq",
        "--fastq2", path + "test_mate2.fastq",
        "-ow"
    ]
    result = runner.invoke(run_docker, args, prog_name='dreem-test')
          
    assert os.path.isfile("output/BitVector_Files/mttr-6-alt-h3_1_134_pop_avg.html")
    util.safe_rmdir("output")
    util.safe_rmdir("input")

    os.chdir('..')
    util.safe_rmdir( DIR_WITH_SPACES )

