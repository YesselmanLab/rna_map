"""
testing the parameters module
"""

import pytest
import os

from rna_map.parameters import parse_parameters_from_file

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
BASE_DIR = os.path.dirname(TEST_DIR)


def test_parse_parameters_empty_file():
    """
    test parse_parameters_from_file with an empty file
    should fill in with all defaults
    """
    path = TEST_DIR + "/resources/test.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == False


def test_parse_parameter_from_default():
    """
    testing the default parameters
    """
    path = TEST_DIR + "/resources/default.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == False


def test_parse_parameter_overwrite():
    path = TEST_DIR + "/resources/changed_params.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == False
    # parameters that were changed
    assert params["bit_vector"]["num_of_surbases"] == 5


def test_partially_filled_parameters():
    path = TEST_DIR + "/resources/few.yml"
    params = parse_parameters_from_file(path)
    assert params["map"]["skip"] == True


def test_invalid_parameters():
    path = TEST_DIR + "/resources/invalid.yml"
    with pytest.raises(Exception):
        parse_parameters_from_file(path)
