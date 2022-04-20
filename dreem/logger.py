import logging
import sys


APP_LOGGER_NAME = 'dreem'

def setup_applevel_logger(logger_name=APP_LOGGER_NAME,
                          is_debug=True,
                          file_name=None):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger

def get_logger(module_name):
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


