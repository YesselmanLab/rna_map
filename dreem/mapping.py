from dreem.parameters import *
from dreem.logger import *
from dreem.util import *

log = init_logger("mapper.py")


class Mapper(object):
    def __init__(self):
        pass

    def run(self, p: Parameters):
        self._p = p
        # run fastqc
        if self._p.map.skip_fastqc:
            self.__skip_method_by_user("fastqc", "skip_fastqc")
            return
        self.__run_fastqc()

    def __skip_method_by_user(self, method_name, method_param):
        log.info(
            "SKIPPING {}, was requested by user using param {}".format(
                method_name, method_param
            )
        )

    def __skip_without_overwrite(self, method_name):
        log.info(
            "SKIPPING {}, it has been run already! specify -overwrite to rerun".format(
                method_name
            )
        )

    # TODO write log to file
    def __run_command(self, method_name, cmd):
        log.info("running {}".format(method_name))
        log.debug(cmd)
        output, error_msg = run_command(cmd)
        if error_msg is not None:
            log_error_and_exit(log, "fastqc", error)
        else:
            log.info("{} ran without errors".format(method_name))

    def __run_fastqc(self):
        fastqc_dir = self._p.dirs.mapping + self._p.ins.fastq1_name + "_fastqc"
        if os.path.isdir(fastqc_dir) and not self._p.overwrite:
            self.__skip_without_overwrite("fastqc")
            return
        if self._p.paired:
            fastqc_cmd = "fastqc --extract {fq1} {fq2} --outdir={dir}".format(
                dir=self._p.dirs.mapping,
                fq1=self._p.ins.fastq1,
                fq2=self._p.ins.fastq2,
            )
        else:
            fastqc_cmd = "fastqc --extract {fq1} -outdir={dir}".format(
                dir=self._p.dirs.mapping, fq1=self._p.ins.fastq1
            )
        self.__run_command("fastqc", fastqc_cmd)
