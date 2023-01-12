import pytest
import yaml
import os

from dreem.parameters import Inputs, parse_parameters_from_file

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
BASE_DIR = os.path.dirname(TEST_DIR)


def test_parse_parameters_empty_file():
    path = TEST_DIR + "/resources/test.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == False

def test_parse_parameter_from_default():
    path = TEST_DIR + "/resources/default.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == False

def test_parse_parameter_overwrite():
    path = TEST_DIR + "/resources/changed_params.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == False
    # parameters that were changed
    assert params["bit_vector"]["qscore_cutoff"] == 15
    assert params["bit_vector"]["num_of_surbases"] == 5
