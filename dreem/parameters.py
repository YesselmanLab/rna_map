import os
import logging
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
            self.ref_fasta = paths[0]  # fasta input file
            self.fastq1 = paths[1]  # fastq1 input file
            self.fastq2 = paths[2]  # fastq2 input file
            self.csv = None

        def __get_name(self, path):
            fname = path.split("/")[-1]
            fname = ".".join(fname.split(".")[:-1])
            return fname

        @property
        def ref_fasta_name(self):
            return self.__get_name(self.ref_fasta)

        @property
        def fastq1_name(self):
            return self.__get_name(self.fastq1)

        @property
        def fastq2_name(self):
            return self.__get_name(self.fastq2)

    class _Dirs(object):
        def __init__(self):
            self.resources = settings.get_py_path() + "/resources/"
            self.input = "input/"
            self.output = "output/"
            self.log = "log/"
            self.mapping = self.output + "Mapping_Files/"
            self.bitvector = self.output + "BitVector_Files/"
            self.cluster = self.output + "/"

    class _Files(object):
        def __init__(self, dirs, inputs):
            self.qhred_ascii = dirs.resources + "/phred_ascii.txt"
            self.tg_fastq1 = dirs.mapping + inputs.fastq1_name + "_val_1.fq"
            self.tg_fastq2 = dirs.mapping + inputs.fastq2_name + "_val_2.fq"
            self.bt2_index = dirs.input + inputs.ref_fasta_name
            self.bt2_alignment_output = dirs.mapping + "aligned.sam"
            self.picard_bam_output = dirs.mapping + "aligned.bam"
            self.picard_sam_output = dirs.mapping + "converted.sam"
            self.picard_sorted_bam_output = dirs.mapping + "aligned_sorted.bam"
            self.picard_metrics_output = dirs.mapping + "metrics.txt"

    class _Map(object):
        def __init__(self):
            self.skip = False
            self.skip_fastqc = False
            self.skip_trim_galore = False
            self.tg_q_cutoff = 20
            self.bt2_alignment_args = (
                "--local,--no-unal,--no-discordant,--no-mixed,-X 1000,-L 12,-p 16"
            )

    class _BitVector(object):
        def __init__(self):
            self.skip = False
            self.qscore_cutoff = 25
            self.num_of_surbases = 10
            self.map_score_cutoff = 15
            self.mutation_count_cutoff = 10
            self.percent_length_cutoff = 0.01
            self.miss_info = "."
            self.ambig_info = "?"
            self.nomut_bit = "0"
            self.del_bit = "1"

    class _EMCluster(object):
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
        if args["csv"] is not None:
            inputs.csv = args["csv"]
        dirs = ParametersFactory._Dirs()
        files = ParametersFactory._Files(dirs, inputs)
        p = Parameters(inputs, dirs, files)
        if args['fastq2'] is not None:
            p.paired = True
        if args['overwrite']:
            log.info(
                "-o/--overwrite supplied, will overwrite previous results with same name"
            )
            p.overwrite = True
        if args['param_file'] is not None:
            self.__parse_param_file(args['param_file'], p)
        if args['log_level'] is not None:
            p.log_level = str_to_log_level(args['log_level'])
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
        self.log_level = logging.INFO
        # parameter groups
        self.map = ParametersFactory._Map()
        self.bit_vector = ParametersFactory._BitVector()
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
    if not os.path.isfile(args["fasta"]):
        raise ValueError(f"fasta file: does not exist {args['fasta']}!")
    else:
        log.info(f"fasta file: {args['fasta']} exists")
    if not os.path.isfile(args['fastq1']):
        raise ValueError("fastq1 file: does not exist {}!".format(args['fastq1']))
    else:
        log.info("fastq file: {} exists".format(args['fastq1']))
    if args['fastq2']:
        if not os.path.isfile(args['fastq2']):
            raise ValueError("fastq2 file: does not exist {}!".format(args['fastq2']))
        else:
            log.info("fastq2 file: {} exists".format(args['fastq2']))
            log.info("two fastq files supplied, thus assuming paired reads")
        return [args['fasta'], args['fastq1'], args['fastq2']]
    else:
        return [args['fasta'], args['fastq1']]


def setup_parameters(args):
    pf = ParametersFactory()
    global parameters
    parameters = pf.get_parameters(args)

    # p.update("map.fastqc", False)
