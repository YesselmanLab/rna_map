import yaml
import json
import jsonschema
from jsonschema import Draft4Validator, validators

from pathlib import Path
from dataclasses import dataclass

from rna_map import settings, logger
from rna_map.settings import get_py_path

log = logger.get_logger("PARAMETERS")


@dataclass(frozen=True, order=True)
class Inputs:
    """
    Input parameters
    """

    fasta: str
    fastq1: str
    fastq2: str = ""
    csv: str = ""

    def is_paired(self):
        """
        Check if the input is paired i.e. has both R1 and R2
        """
        if self.fastq2 != "":
            return True
        return False

    def supplied_csv(self):
        """
        Check if the user supplied a csv file
        """
        if self.csv != "":
            return True
        return False

    def fastq1_name(self):
        """
        Get the name of the fastq1 file
        """
        return Path(self.fastq1).stem

    def fastq2_name(self):
        """
        Get the name of the fastq2 file
        """
        return Path(self.fastq2).stem


def extend_with_default(validator_class):
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property_, subschema in properties.items():
            if "default" in subschema and not isinstance(instance, list):
                instance.setdefault(property_, subschema["default"])

        for error in validate_properties(
            validator,
            properties,
            instance,
            schema,
        ):
            yield error

    return validators.extend(
        validator_class,
        {"properties": set_defaults},
    )


def validate_parameters(params):
    path = get_py_path() + "/resources/params_schema.json"
    with open(path) as f:
        schema = json.load(f)
    # Validate the params against the schema
    FillDefaultValidatingDraft4Validator = extend_with_default(Draft4Validator)
    try:
        FillDefaultValidatingDraft4Validator(schema).validate(params)
    except jsonschema.exceptions.ValidationError as e:
        raise ValueError(e.message)


def parse_parameters_from_file(param_file):
    """
    Parse a YAML file and validate from a schema file loaded from json
    """
    # load param_file and validate against schema
    with open(param_file) as f:
        params = yaml.safe_load(f)
    if params is None:
        params = {}
    validate_parameters(params)
    return params


def get_default_params():
    """
    Get the default parameters
    """
    path = get_py_path() + "/resources/default.yml"
    return parse_parameters_from_file(path)


# abbreviations
# tg -> trim galore
# bt2 -> bowtie 2
