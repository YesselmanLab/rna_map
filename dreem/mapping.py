"""
module to run commands to map sequencing reads to reference sequences
"""

import os

from pathlib import Path
import pandas as pd
from dreem import logger
from dreem.parameters import Inputs
from dreem.settings import get_py_path
from dreem.exception import (
    DREEMInputException,
    DREEMMissingRequirementsException,
    DREEMExternalProgramException,
)
from dreem.external_cmd import (
    does_program_exist,
    run_command,
    get_fastqc_version,
    get_bowtie2_version,
    get_cutadapt_version,
    get_trim_galore_version,
)

log = logger.get_logger("MAPPING")


def skip_without_overwrite(method_name):
    """
    Skip a method if the output directory already exists and overwrite is False
    """
    log.info(
        f"SKIPPING {method_name}, it has been run already! specify -overwrite "
        f"to rerun"
    )


def skip_method_by_user(method_name, method_param):
    """
    Skip a method if the user requests it
    """
    log.info(
        f"SKIPPING {method_name}, was requested by user using param "
        f"{method_param}"
    )



class Mapper(object):
    def check_program_versions(self):
        """
        check to make sure all required programs exist and have the
        right versions
        """
        programs = "fastqc,bowtie2,trim_galore,cutadapt".split(",")
        for prog in programs:
            if not does_program_exist(prog):
                self.__program_not_found(prog)

        log.info(f"bowtie2 {get_bowtie2_version()} detected!")
        log.info(f"fastqc {get_fastqc_version()} detected!")
        log.info(f"trim_galore {get_trim_galore_version()} detected!")
        log.info(f"cutapt {get_cutadapt_version()} detected!")

    def run(self, ins: Inputs, params):
        pass
        # don't rerun unless asked with -overwrite
        # if os.path.isfile(self._p.files.tg_fastq1) and not self._p.map.overwrite:
        #    self.__skip_without_overwrite("trim_galore")
        #    return
        #run_fastqc(ins.fastq1, ins.fastq2, out_dir + "/Mapping_Files/")

        # skip by user request, make sure to copy the files for next step
        # skip trim galore
        """
        if skip:
            skip_method_by_user("trim_glore", "skip_trim_galore")
            log.info("copying input fastq files directly")
            shutil.copy(self._p.ins.fastq1, self._p.files.tg_fastq1)
            if self._p.paired:
                shutil.copy(self._p.ins.fastq2, self._p.files.tg_fastq2)
            return
        """
        # self.__run_trim_glore()  # run trim galore
        # self.__run_bowtie_build()  # run bowtie build
        # self.__run_bowtie_alignment()  # row bowtie
        # self.__run_picard_bam_convert()  # convert sam to bam
        # self.__run_picard_sort()
        # self.__run_picard_metrics()
        log.info("finished mapping!")

    def __program_not_found(self, p_name):
        msg = f"{p_name} is not found make sure it is accessible in $PATH"
        raise DREEMMissingRequirementsException(msg)

    def __skip_method_by_user(self, method_name, method_param):
        log.info(
            "SKIPPING {}, was requested by user using param {}".format(
                method_name, method_param
            )
        )

    def __log_output(self):
        # f = open(f"{self._p.dirs.log}/{method_name}.log", "w")
        # f.write(output)
        # f.close()
        pass

    # run programs #############################################################
