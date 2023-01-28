import yaml
import cloup
from cloup import option_group, option
from rna_map.logger import get_logger

log = get_logger('CLI_OPTS')

def main_options():
    return option_group(
        "Main arguments",
        "These are the main arguments for the command line interface",
        option(
            "-fa",
            "--fasta",
            type=cloup.Path(exists=True),
            required=True,
            help="The fasta file containing the reference sequences",
        ),
        option(
            "-fq1",
            "--fastq1",
            type=cloup.Path(exists=True),
            required=True,
            help="The fastq file containing the single end reads or the first pair of paired end reads",
        ),
        option(
            "-fq2",
            "--fastq2",
            type=str,
            default="",
            help="The fastq file containing the second pair of paired end reads",
        ),
        option(
            "--dot-bracket",
            type=str,
            default="",
            help="The directory containing the input files",
        ),
        option(
            "-pf",
            "--param-file",
            type=str,
            default=None,
            help="A yml formatted file to specify parameters, see rna_map/resources/default.yml for an example",
        ),
        option(
            "-pp",
            "--param-preset",
            type=str,
            default=None,
            help="run a set of parameters for specific uses like 'barcoded-libraries'",
        ),
    )

def docker_options():
    return option_group(
        "Docker options",
        "These are the options for running the command line interface in a docker container",
        option(
            "--docker",
            is_flag=True,
            help="Run the program in a docker container",
        ),
        option(
            "--docker-image",
            type=str,
            default="rna-map",
            help="The docker image to use",
        ),
        option(
            "--docker-platform",
            type=str,
            default="",
            help="The platform to use for the docker image",
        ),
    )

def mapping_options():
    return option_group(
        "Mapping options",
        "These are the options for pre processing of fastq files and alignment to reference sequences",
        option(
            "--skip-fastqc",
            is_flag=True,
            help="do not run fastqc for quality control of sequence data",
        ),
        option(
            "--skip-trim-galore",
            is_flag=True,
            help="do not run trim galore for quality control of sequence data",
        ),
        option(
            "--tg-q-cutoff",
            type=int,
            default=20,
            help="the quality cutoff for trim galore",
        ),
        option(
            "--bt2-alignment-args",
            help="the arguments to pass to bowtie2 for alignment seperated by commas",
        ),
        option(
            "--save-unaligned",
            is_flag=True,
            help="the path to save unaligned reads to",
        ),
    )

def bit_vector_options():
    return option_group(
        "Bit vector options",
        "These are the options for the bit vector step",
        option(
            "--skip-bit-vector",
            is_flag=True,
            help="do not run the bit vector step",
        ),
        option(
            "--summary-output-only",
            is_flag=True,
            help="do not generate bit vector files or plots recommended when there are thousands of reference sequences",
        ),
        option(
            "--plot-sequence",
            is_flag=True,
            help="plot sequence and structure is supplied under the population average plots",
        ),
        option(
            "--map-score-cutoff",
            type=int,
            default=15,
            help="reject any bit vector where the mapping score for bowtie2 alignment is less than this value",
        ),
        option(
            "--qscore-cutoff",
            type=int,
            default=25,
            help="quality score of read nucleotide, sets to ambigious if under this val",
        ),
        option(
            "--mutation-count-cutoff",
            type=int,
            default=5,
            help="maximum number of mutations allowed in a bit vector will be discarded if higher",
        ),
        option(
            "--percent-length-cutoff",
            type=float,
            default=0.1,
            help="minium percent of the length of the reference sequence allowed in a bit vector will be discarded if lower",
        ),
        option(
            "--min-mut-distance",
            type=int,
            default=5,
            help="minimum distance between mutations in a bit vector will be discarded if lower",
        ),
    )

def misc_options():
    return option_group(
        "Misc options",
        "These are the options for the misc stage",
        option(
            "--overwrite",
            is_flag=True,
            help="overwrite the output directory if it exists",
        ),
        option(
            "--restore-org-behavior",
            is_flag=True,
            help="restore the original behavior of the rna_map",
        ),
        option(
            "--stricter-bv-constraints",
            is_flag=True,
            help="use stricter constraints for bit vector generation, use at your own risk!",
        ),
        option(
            "--debug",
            is_flag=True,
            help="enable debug mode",
        ),
    )



def parse_cli_args(params, args):
    # main options
    # docker options
    # mapping options
    if args['skip_fastqc']:
        log.info("skipping fastqc for quality control only do this if you are confident in the quality of your data")
        params["map"]["skip_fastqc"] = args['skip_fastqc']
    if args['skip_trim_galore']:
        log.info("skipping trim galore for quality control not recommended")
        params["map"]["skip_trim_galore"] = args['skip_trim_galore']
    if args['tg_q_cutoff'] != 20:
        log.info("trim galore quality cutoff set to {value}".format(value=args['tg_q_cutoff']))
        params["map"]["tg_q_cutoff"] = args['tg_q_cutoff']
    if args['bt2_alignment_args'] is not None:
        log.info("bowtie2 alignment arguments set to {value}".format(value=args['bt2_alignment_args']))
        params["map"]["bt2_alignment_args"] = args['bt2_alignment_args']
    if args['save_unaligned']:
        log.info("saving unaligned reads to {value}".format(value=args['save_unaligned']))
        params["map"]["save_unaligned"] = args['save_unaligned']
    # bit_vector options
    if args['skip_bit_vector']:
        log.info("skipping bit vector step")
        params["bit_vector"]["skip"] = args['skip_bit_vector']
    if args['summary_output_only']:
        log.info("only outputting summary files")
        params["bit_vector"]["summary_output_only"] = args['summary_output_only']
    if args['plot_sequence']:
        log.info("plotting sequence/structure on bit vector plots")
        params["bit_vector"]["plot_sequence"] = args['plot_sequence']
    if args['map_score_cutoff'] != 15:
        log.info("mapping score cutoff set to {value}".format(value=args['map_score_cutoff']))
        params["bit_vector"]["map_score_cutoff"] = args['map_score_cutoff']
    if args['qscore_cutoff'] != 25:
        log.info("qscore cutoff set to {value}".format(value=args['qscore_cutoff']))
        params["bit_vector"]["qscore_cutoff"] = args['qscore_cutoff']
    if args['mutation_count_cutoff'] != 5:
        log.info("mutation count cutoff set to {value} this will only run if --stricter-bv-constraints is set".format(value=args['mutation_count_cutoff']))
        params["bit_vector"]["stricter_constraints"]["mutation_count_cutoff"] = args['mutation_count_cutoff']
    if args['percent_length_cutoff'] != 0.1:
        log.info("percent length cutoff set to {value} this will only run if --stricter-bv-constraints is set".format(value=args['percent_length_cutoff']))
        params["bit_vector"]["stricter_constraints"]["percent_length_cutoff"] = args['percent_length_cutoff']
    if args['min_mut_distance'] != 5:
        log.info("minimum mutation distance set to {value} this will only run if --stricter-bv-constraints is set".format(value=args['min_mut_distance']))
        params["bit_vector"]["stricter_constraints"]["min_mut_distance"] = args['min_mut_distance']
    # misc options
    if args['overwrite']:
        log.info("will overwrite all existing files")
        params["overwrite"] = args['overwrite']
    if args['restore_org_behavior']:
        log.info("restoring original behavior of rna_map publications")
        params["restore_org_behavior"] = args['restore_org_behavior']
    if args['stricter_bv_constraints']:
        log.info("stricter bit vector constraints are active please use at your own risk")
        params["stricter_bv_constraints"] = args['stricter_bv_constraints']
