import pytest
import yaml
import os
from click.testing import CliRunner
from dreem import parameters, settings, run, run_docker
import dreem

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
BASE_DIR = os.path.dirname(TEST_DIR)


def get_default_args():
    # TODO joe I fixed these params so it would run but idk if these
    # are correct -CJ
    p = {
        'fasta'                : TEST_DIR + '/resources/case_1/test.fasta',
        'fastq1'               : TEST_DIR + '/resources/case_1/test_mate1.fastq',
        'fastq2'               : TEST_DIR + '/resources/case_1/test_mate2.fastq',
        'overwrite'            : False,
        'csv'                  : None,
        'param_file'           : None,
        'log_level'            : 'INFO',
        'bv_overwrite'         : False,
        'restore_org_behavior' : False,
        'bt2_alignment_args'   : None,
        'qscore_cutoff'        : None,
        'num_of_surbases'      : None,
        'map_score_cutoff'     : None,
        'percent_length_cutoff': None,
        'mutation_count_cutoff': None,
        'skip_fastqc'          : False,
        'skip_trim_galore'     : False,
        'plot_sequence'        : False,
        'summary_output_only'  : False
    }
    return p


def test_generate_parameters():
    p = get_default_args()
    params = parameters.ParametersFactory().get_parameters(p)
    assert p['fasta'] == params.ins.ref_fasta


def test_parse_param_file():
    p = get_default_args()
    p['param_file'] = BASE_DIR + "/dreem/resources/default.yml"
    params = parameters.ParametersFactory().get_parameters(p)
    assert params.map.skip == False
    assert params.bit_vector.mutation_count_cutoff == 10


def test_parse_new_param_file():
    p = get_default_args()
    p['param_file'] = TEST_DIR + "/resources/test.yml"
    params = parameters.ParametersFactory().get_parameters(p)
    assert params.map.bt2_alignment_args == \
           "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 20,-p 16"


def test_write_to_yaml():
    p = get_default_args()
    params = parameters.ParametersFactory().get_parameters(p)
    params.to_yaml_file("test.yml")
    loaded_p = yaml.load(open("test.yml"), Loader=yaml.FullLoader)
    k, v = next(iter(loaded_p['bit_vector'][0].items()))
    assert k == "skip"
    os.remove("test.yml")


def test_get_sub_params():
    pf = parameters.ParametersFactory()
    map = pf._Map()
    assert map.skip == False


@pytest.mark.integration
def test_help_strings():
    open('loc', 'w').write(dreem.__file__)
    runner = CliRunner()

    args = [
        "--help"
    ]
    result1 = runner.invoke(run, args, prog_name='dreem-test')
    result2 = runner.invoke(run_docker, args, prog_name='dreem-test')
    assert len(result1.output) > 0
    assert result1.output == result2.output
