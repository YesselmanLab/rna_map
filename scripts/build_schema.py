"""
script for building a schema for yml parameter file
"""
import yaml
import json
import click


def assign_properties(name, schema):
    pass


@click.command()
@click.argument("param_file", type=click.Path(exists=True))
@click.option("--defaults")
@click.option("--required")
def main(param_file, defaults, required):
    """
    main function for script
    """
    with open(param_file, "r") as f:
        params = yaml.safe_load(f)
    assign_properties(params)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
