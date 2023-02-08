"""
Handles logging for module
"""

import logging
import sys

APP_LOGGER_NAME = "rna_map"


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=False, file_name=None):
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

