"""
command line iterface for rna_map
"""
import os
import sys
import subprocess
import cloup

from rna_map.cli_opts import (
    bit_vector_options,
    docker_options,
    main_options,
    mapping_options,
    misc_options,
    parse_cli_args,
)
from rna_map.external_cmd import does_program_exist
from rna_map.logger import setup_applevel_logger, get_logger
from rna_map.parameters import parse_parameters_from_file, get_default_params
from rna_map.run import run
from rna_map.settings import get_py_path

log = get_logger("CLI")


def get_logo() -> str:
    with open(os.path.join(get_py_path(), "resources", "logo.txt"), "r") as f:
        return "".join(f.readlines())


# docker ######################################################################


def check_docker_image(name: str) -> bool:
    try:
        subprocess.check_output(f"docker image inspect {name}", shell=True)
    except subprocess.CalledProcessError:
        return False
    return True


def run_in_docker(args):
    if not does_program_exist("docker"):
        raise ValueError("docker is not installed")
    if not check_docker_image(args["docker_image"]):
        raise ValueError(f"{args['docker_image']} docker image not found")

    file_map = {
        "dot_bracket": "test.csv",
        "param_file": "test.yml",
        "fasta": "test.fasta",
        "fastq1": "test_mate1.fastq",
        "fastq2": "test_mate2.fastq",
    }
    file_args = {
        "dot_bracket": "--dot-bracket",
        "param_file": "--param-file",
        "fasta": "--fasta",
        "fastq1": "--fastq1",
        "fastq2": "--fastq2",
    }
    files = "fasta,fastq1,fastq2,dot_bracket,param_file".split(",")
    platform = args["docker_platform"]
    platform_str = ""
    if platform != "":
        platform_str = "--platform " + platform
    docker_cmd = f"docker run --name rna-map-cont " f"{platform_str} -v $(pwd):/data "
    dreem_cmd = "rna-map "
    dirs = {os.getcwd(): "/data"}
    dcount = 2
    pos = 1
    skip_dreem_args = [
        "-fa",
        "--fasta",
        "-fq1",
        "--fastq1",
        "-fq2",
        "--fastq2",
        "-db",
        "--dot-bracket",
        "-pf",
        "--param-file",
        "--docker-image",
        "--docker-platform",
    ]
    keep_args = []
    while pos < len(sys.argv):
        if sys.argv[pos] in skip_dreem_args:
            pos += 2
            continue
        if sys.argv[pos] == "--docker":
            pos += 1
            continue
        keep_args.append(sys.argv[pos])
        pos += 1
    dreem_cmd += " ".join(keep_args) + " "
    for f in files:
        f_path = args[f]
        if f_path is None or f_path == "":
            continue
        dir_name = os.path.abspath(os.path.dirname(f_path))
        if dir_name == os.path.abspath(os.getcwd()):
            dreem_cmd += f"{file_args[f]} {file_map[f]} "
            continue
        if dir_name not in dirs:
            dirs[dir_name] = f"/data{dcount}"
            docker_cmd += f"-v {dir_name}:/{dirs[dir_name]} "
            dcount += 1
        dreem_cmd += f"{file_args[f]} {dirs[dir_name]}/{file_map[f]} "
    docker_cmd += f"{args['docker_image']} "
    docker_cmd += dreem_cmd
    log.info(
        f" RUNNING DOCKER #####################################################\n"
        f": {docker_cmd}"
    )
    subprocess.call(docker_cmd, shell=True)
    subprocess.call("docker rm rna-map-cont", shell=True)


# cli #########################################################################


@cloup.command()
@main_options()
@mapping_options()
@bit_vector_options()
@docker_options()
@misc_options()
def cli(**args):
    """
    rapid analysis of RNA mutational profiling (MaP) experiments.
    """
    # setup logging
    if args["debug"]:
        setup_applevel_logger(is_debug=True)
        log.info("Debug logging is on")
    else:
        setup_applevel_logger()
    # check to see if we are running in docker
    if args["docker"]:
        print("running in docker")
        run_in_docker(args)
        return
    # log commandline options
    log.info("\n" + get_logo())
    log.info("ran at commandline as: ")
    log.info(" ".join(sys.argv))
    # setup parameters
    if args["param_file"] is not None and args["param_preset"] is not None:
        raise ValueError("cannot specify both param_file and param_preset")
    if args["param_preset"] is not None:
        param_path = os.path.join(
            get_py_path(), "resources", "presets", f"{args['param_preset']}.yml"
        )
        if not os.path.isfile(param_path):
            raise ValueError(f"preset {args['param_preset']} does not exist")
        log.info(f"using param preset: {args['param_preset']}")
        params = parse_parameters_from_file(param_path)
    elif args["param_file"] is not None:
        params = parse_parameters_from_file(args["param_file"])
    else:
        params = get_default_params()
    parse_cli_args(params, args)
    # run rna_map pipeline
    run(
        args.pop("fasta"),
        args.pop("fastq1"),
        args.pop("fastq2"),
        args.pop("dot_bracket"),
        params,
    )


if __name__ == "__main__":
    cli()
