import os
import logging
import colorlog

class DreemError(Exception):
    pass

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
                "DEBUG": "cyan",
                "WARNING": "yellow",
                "ERROR": "red",
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


def log_error_and_exit(log, msg):
    log.error(msg)
    raise DreemError(msg)


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

log = init_logger("bit_vector.py", "dreem.log", start=True)
