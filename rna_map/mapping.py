"""
module to run commands to map sequencing reads to reference sequences
"""

import os
import shutil
from typing import Dict
from pathlib import Path
from rna_map import logger
from rna_map.parameters import Inputs
from rna_map.exception import (
    DREEMMissingRequirementsException,
)
from rna_map.external_cmd import (
    does_program_exist,
    get_fastqc_version,
    get_bowtie2_version,
    get_cutadapt_version,
    get_trim_galore_version,
    run_bowtie_alignment,
    run_bowtie_build,
    run_fastqc,
    run_trim_glore,
)

log = logger.get_logger("MAPPING")


class Mapper(object):
    def __init__(self):
        self.__setup = False

    def setup(self, params: Dict):
        self.__params = params
        self.__setup = True
        # params
        self.__in_dir = self.__params["dirs"]["input"]
        self.__out_dir = os.path.join(
            self.__params["dirs"]["output"], "Mapping_Files"
        )
        self.__overwrite = self.__params["overwrite"]
        self.check_program_versions()
        self.build_directories(self.__params)

    def build_directories(self, params):
        log.info("building directory structure")
        os.makedirs(params["dirs"]["log"], exist_ok=True)
        os.makedirs(params["dirs"]["input"], exist_ok=True)
        os.makedirs(params["dirs"]["output"], exist_ok=True)
        os.makedirs(
            os.path.join(params["dirs"]["output"], "Mapping_Files"),
            exist_ok=True,
        )

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

    def run(self, ins: Inputs):
        if not self.__setup:
            raise Exception("Must call setup before running")
        self.__ins = ins
        self.__run_fastqc()
        self.__run_trim_glore()
        self.__run_bowtie_build()  # run bowtie build
        self.__run_bowtie_alignment()  # row bowtie
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

    def skip_without_overwrite(self, method_name):
        """
        Skip a method if the output directory already exists and overwrite is False
        """
        log.info(
            f"SKIPPING {method_name}, it has been run already! specify -overwrite "
            f"to rerun"
        )

    def skip_method_by_user(self, method_name, method_param):
        """
        Skip a method if the user requests it
        """
        log.info(
            f"SKIPPING {method_name}, was requested by user using param "
            f"{method_param}"
        )

    # run programs #############################################################
    def __run_fastqc(self):
        """
        run fastqc on the fastq files
        """
        if self.__params["map"]["skip_fastqc"]:
            self.skip_method_by_user("fastqc", "skip_fastqc")
            return
        # don't rerun unless asked with -overwrite
        if (
            os.path.isdir(os.path.join(self.__out_dir, "fastqc"))
            and not self.__overwrite
        ):
            self.skip_without_overwrite("fastqc")
            return
        run_fastqc(self.__ins.fastq1, self.__ins.fastq2, self.__out_dir)

    def __run_trim_glore(self):
        """
        run trim galore on the fastq files
        """
        # version
        fq1_path = f"{self.__out_dir}/{Path(self.__ins.fastq1).stem}_val_1.fq"
        if self.__params["map"]["skip_trim_galore"]:
            self.skip_method_by_user("trim_glore", "skip_trim_galore")
            shutil.copy(self.__ins.fastq1, fq1_path)
            fq2_path = (
                f"{self.__out_dir}/{Path(self.__ins.fastq2).stem}_val_2.fq"
            )
            shutil.copy(self.__ins.fastq2, fq2_path)
            return
        # don't rerun unless asked with -overwrite
        if os.path.isfile(fq1_path) and not self.__overwrite:
            self.skip_without_overwrite("trim_glore")
            return
        return run_trim_glore(
            self.__ins.fastq1, self.__ins.fastq2, self.__out_dir
        )

    def __run_bowtie_build(self):
        """
        run bowtie build on the reference sequence
        """
        bt2_index = os.path.join(
            self.__in_dir, f"{Path(self.__ins.fasta).stem}.bt2"
        )
        if os.path.isfile(bt2_index) and not self.__overwrite:
            self.skip_without_overwrite("bowtie_build")
            return
        return run_bowtie_build(self.__ins.fasta, self.__in_dir)

    def __run_bowtie_alignment(self):
        """
        run bowtie alignment on the fastq files
        """
        sam_path = os.path.join(self.__out_dir, "aligned.sam")
        if os.path.isfile(sam_path) and not self.__overwrite:
            self.skip_without_overwrite("bowtie_alignment")
            return
        if self.__ins.is_paired():
            fq1_path = (
                f"{self.__out_dir}/{Path(self.__ins.fastq1).stem}_val_1.fq"
            )
            fq2_path = (
                f"{self.__out_dir}/{Path(self.__ins.fastq2).stem}_val_2.fq"
            )
        else:
            fq1_path = (
                f"{self.__out_dir}/{Path(self.__ins.fastq1).stem}_trimmed.fq"
            )
            fq2_path = ""
        return run_bowtie_alignment(
            self.__ins.fasta,
            fq1_path,
            fq2_path,
            self.__in_dir,
            self.__out_dir,
            self.__params["map"]["bt2_alignment_args"],
        )