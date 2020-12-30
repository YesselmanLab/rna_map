import os
import subprocess


def safe_mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


def run_command(cmd):
    output, error_msg = None, None
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT
        ).decode("utf8")
    except subprocess.CalledProcessError as exc:
        error_msg = exc.output.decode("utf8")
    return output, error_msg


def log_error_and_exit(log, pname, error_msg):
    log.error("{} returned error:".format(pname))
    log.error(error_msg)
    log.error("EXITING")
    exit()
