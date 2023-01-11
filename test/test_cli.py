
from dreem.run import parse_parameters_from_file

def test_parse_parameters_from_file():
    path = "test.yml"
    params = parse_parameters_from_file(path)
    print(params)