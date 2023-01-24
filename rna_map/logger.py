"""
Handles logging for module
"""

import logging
import sys

APP_LOGGER_NAME = "rna_map"


def setup_applevel_logger(
    logger_name=APP_LOGGER_NAME, is_debug=False, file_name=None
):
    """
    Set up the logger for the app
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    # pylint: disable=C0103
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(module_name):
    """
    Get the logger for the module
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)


def log_error_and_exit(log, msg, exception):
    log.error(msg)
    raise exception(msg)


def str_to_log_level(s: str):
    s = s.rstrip().lstrip().lower()
    if s == "info":
        return logging.INFO
    elif s == "debug":
        return logging.DEBUG
    elif s == "warn":
        return logging.WARN
    elif s == "error":
        return logging.ERROR
    elif s == "critical":
        return logging.CRITICAL
    else:
        raise ValueError("unknown log level: {}".format(s))
