import os

from dreem import logger

app_logger = logger.setup_applevel_logger()

def test_logging(caplog):
    app_logger.info("test")
    assert(caplog.record_tuples[0][0] == 'dreem')
    assert(caplog.record_tuples[0][-1] == 'test')

def test_logging_2(caplog):
    log = logger.get_logger('TEST')
    log.info("test")
    assert(caplog.record_tuples[0][0] == 'dreem.TEST')
    assert(caplog.record_tuples[0][-1] == 'test')

def test_write_to_file():
    app_logger = logger.setup_applevel_logger(file_name='dreem.log')
    app_logger.info("test")
    assert(os.path.isfile("dreem.log"))
    os.remove("dreem.log")
