import click
import pandas as pd

from dreem import logger, run, run_docker

log = logger.log


def get_default_args():
    p = {
        'fasta'                : None,
        'fastq1'               : None,
        'fastq2'               : None,
        'db'                   : None,
        'param_file'           : None,
        'overwrite'            : False,
        'log_level'            : 'INFO',
        'restore_org_behavior' : False,
        'map_overwrite'        : False,
        'skip'                 : False,
        'skip_fastqc'          : False,
        'skip_trim_galore'     : False,
        'bt2_alignment_args'   : None,
        'bv_overwrite'         : False,
        'qscore_cutoff'        : None,
        'num_of_surbases'      : None,
        'map_score_cutoff'     : None,
        'mutation_count_cutoff': None,
        'percent_length_cutoff': None}
    return p


def set_args_from_row(p, row):
    pass


@click.command()
@click.argument(csv)
def main(csv):
    df = pd.read_csv(csv)
    if 'fasta' not in df:
        logger.log_error_and_exit(log, "'fasta' must be in a row in the csv file")
    for i, row in df.iterrows():
        pass


if __name__ == "__main__":
    main()
