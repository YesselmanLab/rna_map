import os
import shutil

from dreem.parameters import *
from dreem.logger import *
from dreem.util import *

log = init_logger("mapper.py", "dreem.log")

class DreemMappingError(Exception):
    pass

class Mapper(object):
    def __init__(self):
        # check programs
        programs = "fastqc,bowtie2,trim_galore,cutadapt".split(",")
        for prog in programs:
            if not does_program_exist(prog):
                self.__program_not_found(log, prog)

        log.info("bowtie2 {} detected!".format(get_bowtie2_version()))
        log.info("fastqc {} detected!".format(get_fastqc_version()))
        log.info("trim_galore {} detected!".format(get_trim_galore_version()))
        log.info("cutapt {} detected!".format(get_cutadapt_version()))

    def run(self, p: Parameters):
        self._p = p
        log.setLevel(p.log_level)
        self.__run_fastqc()  # run fastqc
        self.__run_trim_glore()  # run trim galore
        self.__run_bowtie_build()  # run bowtie build
        self.__run_bowtie_alignment()  # row bowtie
        self.__run_picard_bam_convert()  # convert sam to bam
        self.__run_picard_sort()
        # self.__run_picard_metrics()
        log.info("finished mapping!")

    def __program_not_found(self, log, p_name):
        log_error_and_exit(
                log, "{} is not found make sure it is accessible in $PATH".format(p_name)
        )

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

    def __run_command(self, method_name, cmd):
        log.info("running {}".format(method_name))
        log.debug(cmd)
        output, error_msg = run_command(cmd)
        if error_msg is not None:
            self.__log_error_msg_and_exit(log, method_name, error_msg)
        else:
            log.info("{} ran without errors".format(method_name))
        f = open(f"{self._p.dirs.log}/{method_name}.log", "w")
        f.write(output)
        f.close()
        return output, error_msg

    def __log_error_msg_and_exit(self, log, pname, error_msg):
        s = "{} returned error:".format(pname)
        s = s + "\n" + error_msg + "\n" + "EXITING"
        logger.log_error_and_exit(log, s)

    # TODO check if correct arg val?
    def __validate_bt2_args(self, args):
        f = open(self._p.dirs.resources + "bowtie2.args")
        lines = f.readlines()
        f.close()
        valid_bt2_args = {}
        for l in lines:
            arg, arg_type = l.rstrip().split(",")
            valid_bt2_args[arg] = arg_type
        supplied_args = self._p.map.bt2_alignment_args.split(";")
        for full_arg in supplied_args:
            if full_arg in valid_bt2_args:
                log.debug("{} is a valid bt2 argument".format(full_arg))
                continue
            spl = full_arg.split()
            if len(spl) == 0:
                log_error_and_exit(
                        log, "{} is an invalid bt2 argument".format(full_arg)
                )
            arg, arg_val = spl[0], spl[1]
            if arg in valid_bt2_args:
                log.debug("{} is a valid bt2 argument".format(full_arg))
            else:
                log_error_and_exit(
                        log, "{} is an invalid bt2 argument".format(full_arg)
                )
        log.debug("all bt2 arguments are valid")

    # run programs
    def __run_fastqc(self):
        if self._p.map.skip_fastqc:
            self.__skip_method_by_user("fastqc", "skip_fastqc")
            return
        fastqc_dir = self._p.dirs.mapping + self._p.ins.fastq1_name + "_fastqc"
        if os.path.isdir(fastqc_dir) and not self._p.map.overwrite:
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

    def __run_trim_glore(self):
        # skip by user request, make sure to copy the files for next step
        if self._p.map.skip_trim_galore:
            self.__skip_method_by_user("trim_glore", "skip_trim_galore")
            log.info("copying input fastq files directly")
            shutil.copy(self._p.ins.fastq1, self._p.files.tg_fastq1)
            if self._p.paired:
                shutil.copy(self._p.ins.fastq2, self._p.files.tg_fastq2)
            return
        # don't rerun unless asked with -overwrite
        if os.path.isfile(self._p.files.tg_fastq1) and not self._p.map.overwrite:
            self.__skip_without_overwrite("trim_galore")
            return
        if self._p.paired:
            trim_galore_cmd = "trim_galore --fastqc --paired {fq1} {fq2} -o {out}".format(
                    out=self._p.dirs.mapping, fq1=self._p.ins.fastq1, fq2=self._p.ins.fastq2
            )
        else:
            trim_galore_cmd = "trim_galore --fastqc {fq1} -o {out}".format(
                    out=self._p.dirs.mapping, fq1=self._p.ins.fastq1
            )
        self.__run_command("trim_galore", trim_galore_cmd)

    def __run_bowtie_build(self):
        # TODO assume I already checked is a valid fasta file format
        if os.path.isfile(self._p.dirs.input + self._p.ins.ref_fasta_name + ".1.bt2") and \
                not self._p.overwrite:
            self.__skip_without_overwrite("bowtie-build")
            return

        bowtie_build_cmd = "bowtie2-build \"{fa}\" input/{fa_name}".format(
                fa=self._p.ins.ref_fasta, fa_name=self._p.ins.ref_fasta_name
        )
        self.__run_command("bowtie2-build", bowtie_build_cmd)

    def __run_bowtie_alignment(self):
        if os.path.isfile(self._p.files.bt2_alignment_output) and not self._p.overwrite:
            self.__skip_without_overwrite("bowtie2 alignment")
            return
        # check to make sure bt2 args are valid
        self.__validate_bt2_args(self._p.map.bt2_alignment_args)
        bt2_args = " ".join(self._p.map.bt2_alignment_args.split(";"))
        if self._p.paired:
            bowtie2_cmd = (
                    "bowtie2 " + bt2_args + " -x {btindex} -1 {fq1} -2 {fq2} -S {samfile}"
            ).format(
                    btindex=self._p.files.bt2_index,
                    fq1=self._p.files.tg_fastq1,
                    fq2=self._p.files.tg_fastq2,
                    samfile=self._p.files.bt2_alignment_output,
            )
        else:
            bowtie2_cmd = (
                    "bowtie2 " + bt2_args + " -x {btindex} -U {fq} -S {samfile}"
            ).format(
                    btindex=self._p.files.bt2_index,
                    fq=self._p.files.tg_fastq1,
                    samfile=self._p.files.bt2_alignment_output,
            )

        output, error_msg = self.__run_command("bowtie2 alignment", bowtie2_cmd)
        output_lines = output.split("\n")
        keep = []
        for l in output_lines:
            if len(l) == 0:
                continue
            if l[0] != "U":
                keep.append(l)
        log.info("results for bowtie alignment: \n" + "\n".join(keep))

    def __run_picard_bam_convert(self):
        if os.path.isfile(self._p.files.picard_bam_output) and not self._p.overwrite:
            self.__skip_without_overwrite("picard BAM conversion")
            return

        log.info("Converting BAM file to SAM file format")
        picard_path = self._p.dirs.resources + "/picard.jar"
        picard_cmd = (
                "java -jar "
                + picard_path
                + " SamFormatConverter I={sf} O={bf}".format(
                sf=self._p.files.bt2_alignment_output,
                bf=self._p.files.picard_bam_output,
        )
        )
        self.__run_command("picard BAM conversion", picard_cmd)

    def __run_picard_sort(self):
        if os.path.isfile(self._p.files.picard_sorted_bam_output) and not self._p.overwrite:
            self.__skip_without_overwrite("picard BAM sort")
            return

        log.info("sorting BAM file")
        picard_path = self._p.dirs.resources + "/picard.jar"
        picard_cmd = (
                "java -jar "
                + picard_path
                + " SortSam I={} O={} SORT_ORDER=coordinate".format(
                self._p.files.picard_bam_output, self._p.files.picard_sorted_bam_output
        )
        )
        self.__run_command("picard BAM sort", picard_cmd)

    def __run_picard_metrics(self):
        # if os.path.isfile(self._p.files.picard_sorted_bam_output):
        #    self.__skip_without_overwrite("picard BAM sort")
        #    return

        log.info("generating metrics on BAM file")
        picard_path = self._p.dirs.resources + "/picard.jar"
        picard_cmd = (
                "java -jar "
                + picard_path
                + " CollectMultipleMetrics I={} O={} R={}".format(
                self._p.files.picard_sorted_bam_output,
                self._p.files.picard_metrics_output,
                self._p.ins.ref_fasta,
        )
        )
        self.__run_command("picard metrics", picard_cmd)
