import click
from click_option_group import optgroup
from dreem import parameters



"""
def run_options(function):
    # required options
    function = click.option("-fa", "--fasta", type=click.Path(exists=True), required=True,
                            help="reference sequences in fasta format")(function)
    function = click.option("-fq1", "--fastq1", type=click.Path(exists=True), required=True,
                            help="fastq sequencing file of mate 1")(function)
    # optional args
    # options
    function = click.option("-fq2", "--fastq2", type=click.Path(exists=True),
                            help="fastq sequencing file of mate 2")(function)
    function = click.option("--db", type=click.Path(exists=True),
                            help="A csv formatted file that contains dot bracket info for "
                                 "each sequence")(function)
    function = click.option("-pf", "--param-file", type=click.Path(exists=True),
                            help="A csv formatted file that contains dot bracket info for "
                                 "each sequence")(function)
    # flags
    function = click.option("-ow", "--overwrite", is_flag=True,
                            help="overwrites previous results, if not set will keep  "
                                 "previous calculation checkpoints")(function)

    pf = parameters.ParametersFactory()
    map_args = pf._Map().__dict__
    for k, v in map_args.items():
        if type(v) == bool:
            function = click.option("--map-"+k, is_flag=True)(function)

    return function


@click.option(
        "-ll", "--log-level", help="set log level", default="INFO"
)
"""
