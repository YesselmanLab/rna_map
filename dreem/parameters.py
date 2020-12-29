import os
from future import types
from dreem import settings
from dreem.logger import *

log = init_logger("parameters.py")

# abbreviations
# tg -> trim galore
# bt2 -> bowtie 2


class ParametersFactory(object):
    def __init__(self):
        pass

    class _Inputs(object):
        def __init__(self, paths):
            self.__ref_fasta = paths[0]  # fasta input file
            self.__fastq1 = paths[1]  # fastq1 input file
            self.__fastq2 = paths[2]  # fastq2 input file

        def __get_name(self, path):
            fname = path.split("/")[-1]
            fname = ".".join(fname.split(".")[:-1])
            return fname

        @property
        def ref_fasta(self):
            return self.__ref_fasta

        @property
        def fastq1(self):
            return self.__fastq1

        @property
        def fastq2(self):
            return self.__fastq2

        @property
        def ref_fasta_name(self):
            return self.__get_name(self.__ref_fasta)

        @property
        def fastq1_name(self):
            return self.__get_name(self.__fastq1)

        @property
        def fastq2_name(self):
            return self.__get_name(self.__fastq2)

    class _Dirs(object):
        def __init__(self):
            self.__resources = settings.get_lib_path() + "/resources/"
            self.__input = "input"
            self.__output = "output"
            self.__log = "log"
            self.__mapping = self.__output + "/Mapping_Files/"
            self.__bitvector = self.__output + "/BitVector_Files/"
            self.__cluster = self.__output + "/"

        @property
        def resources(self):
            return self.__resources

        @property
        def input(self):
            return self.__input

        @property
        def output(self):
            return self.__output

        @property
        def log(self):
            return self.__log

        @property
        def mapping(self):
            return self.__mapping

        @property
        def bitvector(self):
            return self.__bitvector

        @property
        def cluster(self):
            return self.__cluster

    class _Files(object):
        def __init__(self, dirs, inputs):
            self.__qhred_ascii = dirs.resources + "/phred_ascii.txt"
            self.__tg_fastq1 = dirs.mapping + "/" + inputs.fastq1 + "_val_1.fq"
            self.__tg_fastq2 = dirs.mapping + "/" + inputs.fastq2 + "_val_2.fq"
            self.__bt2_index = dirs.input + "/" + inputs.ref_fasta_name
            self.__bt2_alignment_output = dirs.mapping + "aligned.sam"
            self.__picard_bam_output = dirs.mapping + "aligned.bam"
            self.__picard_sam_output = dirs.mapping + "converted.sam"

        @property
        def qhred_ascii(self):
            return self.__qhred_ascii

        @property
        def tg_fastq1(self):
            return self.__tg_fastq1

        @property
        def tg_fastq2(self):
            return self.__tg_fastq2

        @property
        def bt2_index(self):
            return self.__bt2_index

        @property
        def bt2_alignment_output(self):
            return self.__bt2_alignment_output

        @property
        def picard_bam_output(self):
            return self.__picard_bam_output

        @property
        def picard_sam_output(self):
            return self.__picard_sam_output

    class _Map(object):
        def __init__(self):
            self.skip = False
            self.skip_fastqc = False
            self.skip_trim_galore = False
            self.tg_q_cutoff = 20
            self.bt2_alignment_args = (
                "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-p 16"
            )

    class __BitVector(object):
        def __init__(self):
            pass

    def __parse_param_file(self, pf_path, p):
        if not os.path.isfile(pf_path):
            log.error("parameter file {} supplied does not exist!".format(pf_path))
            exit()
        f = open(pf_path)
        lines = f.readlines()
        f.close()
        if len(lines) == 0:
            log.warn("parameter file {} supplied but is empty!".format(pf_path))
        for l in lines:
            spl = [x.rstrip().lstrip() for x in l.split("=")]
            if len(spl) == 1:
                p.update(spl[0], True)
            else:
                p.update(spl[0], spl[1])

    def get_parameters(self, args):
        input_files = validate_files(args)
        inputs = ParametersFactory._Inputs(input_files)
        dirs = ParametersFactory._Dirs()
        files = ParametersFactory._Files(dirs, inputs)
        p = Parameters(inputs, dirs, files)
        if args.fastq2 is not None:
            p.paired = True
        if args.overwrite is not None:
            p.overwrite = True

        if args.params is not None:
            self.__parse_param_file(args.params, p)
        exit()
        return p


class Parameters(object):
    def __init__(
        self,
        inputs: ParametersFactory._Inputs,
        dirs: ParametersFactory._Dirs,
        files: ParametersFactory._Files,
    ):
        # general global options
        self.overwrite = False
        self.paired = False
        # parameter groups
        self.map = ParametersFactory._Map()
        self.ins: ParametersFactory._Inputs = inputs
        self.dirs: ParametersFactory._Dirs = dirs
        self.files: ParametersFactory._Files = files

    def update(self, param, val):
        spl = param.split(".")
        if len(spl) == 1:
            self.__dict__.update({param: val})
        else:
            self.__update_parameter(spl[0], spl[1], val)

    def __update_parameter(self, sub, param, val):
        valid = "map,files,dirs".split(",")
        if sub not in valid:
            log_error_and_exit(log, "unknown parameter group: {}".format(sub))
        sub_obj = self.__dict__[sub]
        if param not in sub_obj.__dict__:
            log_error_and_exit(
                log, "unknown parameter {} in group {}".format(sub, param)
            )
        else:
            sub_obj.__dict__.update({param: val})
            log.info(
                "setting parameter {}.{} to {} specified by parameter file".format(
                    sub, param, val
                )
            )


parameters: Parameters = None


def get_parameters() -> Parameters:
    global parameters
    return parameters


# TODO finish
def validate_fasta_file(fpath):
    pass


def validate_files(args):
    if not os.path.isfile(args.fasta):
        raise ValueError("fasta file: does not exist {}!".format(args.fasta))
    else:
        log.info("fasta file: {} exists".format(args.fasta))
    if not os.path.isfile(args.fastq1):
        raise ValueError("fastq1 file: does not exist {}!".format(args.fastq1))
    else:
        log.info("fastq file: {} exists".format(args.fastq1))
    if args.fastq2:
        if not os.path.isfile(args.fastq2):
            raise ValueError("fastq2 file: does not exist {}!".format(args.fastq2))
        else:
            log.info("fastq2 file: {} exists".format(args.fastq2))
            log.info("two fastq files supplied, thus assuming paired reads")
        return [args.fasta, args.fastq1, args.fastq2]
    else:
        return [args.fasta, args.fastq1]


def setup_parameters(args):
    pf = ParametersFactory()
    global parameters
    parameters = pf.get_parameters(args)

    # p.update("map.fastqc", False)
