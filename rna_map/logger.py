"""
Logging Setup Module

This module provides functions to set up logging configurations for Python applications.
It includes functionalities to set up root-level logging, application-level logging, and
to get logger instances with specific module names.

Functions:
    - setup_logging(file_name: str = None) -> logging.Logger
        Set up the root logging configuration with optional file logging.

    - get_logger(module_name: str = "") -> logging.Logger
        Get a logger instance with the specified module name.
"""

import logging
import sys

# logging #####################################################################

APP_LOGGER_NAME = "rna-map"

def setup_logging(file_name: str = None, debug: bool = False) -> logging.Logger:
    """
    Set up logging configuration.

    Args:
        file_name (str, optional): The name of the log file. If provided, logs will be
        written to this file.

    Returns:
        None
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(
        logging.DEBUG if debug else logging.INFO
    )  # Set the root logger level

    # Create a stream handler for output to console
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(
        logging.DEBUG if debug else logging.INFO
    )  # Set the desired level for console output
    formatter = logging.Formatter("%(levelname)s - %(name)s - %(message)s")
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        root_logger.addHandler(fh)

    return root_logger


def setup_applevel_logger(
    logger_name: str = APP_LOGGER_NAME, is_debug: bool = False, file_name: str = None
) -> logging.Logger:
    """
    Set up and configure an application-level logger.

    Args:
        logger_name (str, optional): The name of the logger. Defaults to APP_LOGGER_NAME.
        is_debug (bool, optional): Flag indicating whether to enable debug level logging.
        Defaults to False.
        file_name (str, optional): The name of the log file. Defaults to None.

    Returns:
        logging.Logger: The configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    # pylint: disable=C0103
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(module_name: str = "") -> logging.Logger:
    """
    Get a logger instance with the specified module name.

    Args:
        module_name (str): The name of the module to be included in the logger name.

    Returns:
        logging.Logger: A logger instance with the specified module name.

    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)


