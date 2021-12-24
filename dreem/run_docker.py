import click
import subprocess
import os
import logging
import colorlog
from click_option_group import optgroup


def init_logger(dunder_name, log_outfile=None, testing_mode=False, start=False) -> logging.Logger:
    log_format = (
        "[%(asctime)s " "%(name)s " "%(funcName)s] " "%(levelname)s " "%(message)s"
    )
    bold_seq = "\033[1m"
    colorlog_format = f"{bold_seq}" "%(log_color)s" f"{log_format}"
    logger = logging.getLogger(dunder_name)
    # colorlog.basicConfig(format=colorlog_format, datefmt="%H:%M")
    handler = colorlog.StreamHandler()
    handler.setFormatter(
            colorlog.ColoredFormatter(
                    colorlog_format,
                    datefmt="%H:%M",
                    reset=True,
                    log_colors={
                        "DEBUG"   : "cyan",
                        "WARNING" : "yellow",
                        "ERROR"   : "red",
                        "CRITICAL": "red,bg_white",
                    },
            )
    )

    logger.addHandler(handler)

    if log_outfile is not None:
        if start:
            if os.path.isfile(log_outfile):
                os.remove(log_outfile)
        fileHandler = logging.FileHandler(log_outfile)
        fileHandler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(fileHandler)

    if testing_mode:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    return logger


log = init_logger('RUN_DOCKER')


def check_docker_image(name):
    try:
        output = subprocess.check_output(f'docker image inspect {name}', shell=True)
    except subprocess.CalledProcessError:
        return False
    return True


@click.command()

@optgroup.group("main arguments")
@optgroup.option("-fa", "--fasta", type=click.Path(exists=True), required=True,
                 help="reference sequences in fasta format")
@optgroup.option("-fq1", "--fastq1", type=click.Path(exists=True), required=True,
                 help="fastq sequencing file of mate 1")
@optgroup.option("-fq2", "--fastq2", type=click.Path(exists=True),
                 help="fastq sequencing file of mate 2", default=None)
@optgroup.group("common options")
@optgroup.option("--dot_bracket", type=click.Path(exists=True),
                help="A csv formatted file that contains dot bracket info for each sequence")
@optgroup.option("-pf", "--param-file", type=click.Path(exists=True),
                help="A yml formatted file to specify parameters")
@optgroup.option("-ow", "--overwrite", is_flag=True,
                help="overwrites previous results, if not set will keep previous " 
                     "calculation checkpoints")
@optgroup.option("-ll", "--log-level", help="set log level (INFO|WARN|DEBUG|ERROR|FATAL)", 
                default="INFO")
@optgroup.option("-rob", "--restore_org_behavior", is_flag=True, default=False, 
                help="retores the original behavior of dreem upon first release")
@optgroup.group("map options")
@optgroup.option("--map-overwrite", is_flag=True,
help="overwrite mapping calculation")
@optgroup.option("--skip", is_flag=True,
help="do not perform sequence mapping, not recommended")
@optgroup.option("--skip_fastqc", is_flag=True,
help="do not run fastqc for quality control of sequence data")
@optgroup.option("--skip_trim_galore", is_flag=True,
help="do not run trim galore to remove adapter sequences at ends")
@optgroup.option("--tg_q_cutoff", default=None,
help="TODO")
@optgroup.option("--bt2_alignment_args", default=None,
help="TODO")
@optgroup.group("bv options")
@optgroup.option("--skip", is_flag=True,
help="skip bit vector generation step, not recommended")
@optgroup.option("--bv-overwrite", is_flag=True,
help="overwrite bit vector calculation")
@optgroup.option("--qscore_cutoff", default=None,
help="quality score of read nucleotide, sets to ambigious if under this val")
@optgroup.option("--num_of_surbases", default=None,
help="number of bases around a mutation")
@optgroup.option("--map_score_cutoff", default=None,
help="map alignment score cutoff for a read, read is discarded if under this value")
@optgroup.option("--mutation_count_cutoff", default=None,
help="maximum number of mutations in a read allowable")
@optgroup.option("--percent_length_cutoff", default=None,
help="read is discarded if less than this percent of a ref sequence is included")
@optgroup.option("--summary_output_only", is_flag=True,
help="")
@optgroup.option("--plot_sequence", is_flag=True,
help="")
def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    # generate docker container and run code
    docker_img = 'dreem'
    if not check_docker_image(docker_img):
        docker_img = 'jyesselm/dreem'
        if not check_docker_image(docker_img):
            log.error("cannot find docker image. Make sure you have it built or downloaded")
            exit()
    cmd = "docker run "
    files = "fasta,fastq1,fastq2,dot_bracket,param_file".split(",")
    file_map = {
        'dot_bracket'        : 'test.csv',
        'param_file': 'test.param',
        'fasta'     : 'test.fasta',
        'fastq1'    : 'test_mate1.fastq',
        'fastq2'    : 'test_mate2.fastq'
    }
    
    abs_files = []
    for f in files:
        f_path = args[f]
        if f_path is None:
            continue
        spl = f_path.split("/")
        abs_files.append(os.path.abspath(f_path) + ":/data/" + file_map[f])
    
    for af in abs_files:
        #cmd += f' -v "{af}" '
        cmd += f' -v {af} '

    cmd += f" --name dreem_cont -t {docker_img}  dreem -fa test.fasta -fq1 test_mate1.fastq -fq2 test_mate2.fastq "
    
    del args['fasta']
    del args['fastq1']
    del args['fastq2']
    
    for k, v in args.items():
        if v is None or v == False:
            continue
        if k in file_map:
            v = file_map[k]
        
        if k.find('dot') == -1:  # added just for dot-bracket
            k = k.replace("_", "-")
        # basically need to check if its a flag or not 
        if not isinstance(v, bool): 
            cmd += f"--{k} '{v}' "
        else:
            cmd += f"--{k} "

    log.info("DOCKER CMD:\n" + cmd)
    subprocess.call(cmd, shell=True)
    # this should probably be wrapped in some kind of a bigger try catch 
    # loop that isolates when something doesn't work
    log.info("clean up and copy files from docker")
    log.info("docker cp dreem_cont:/data/output output")
    log.info("docker cp dreem_cont:/data/log log")
    log.info("docker cp dreem_cont:/data/dreem.log dreem.log")
    log.info("docker rm dreem_cont")
    subprocess.call("docker cp dreem_cont:/data/output output", shell=True)
    subprocess.call("docker cp dreem_cont:/data/log log", shell=True)
    subprocess.call("docker cp dreem_cont:/data/dreem.log dreem.log", shell=True)
    subprocess.call("docker rm dreem_cont", shell=True)
if __name__ == "__main__":
    main()
