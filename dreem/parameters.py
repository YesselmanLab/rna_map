import os
import logging
import yaml
from future import types
from dreem import settings, logger

log = logger.log
log_error_and_exit = logger.log_error_and_exit


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
            if len(paths) == 3:
                self.fastq2 = paths[2]  # fastq2 input file
            else:
                self.fastq2 = None
            self.dot_bracket = None

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
            if inputs.fastq2 is not None:
                self.tg_fastq1 = dirs.mapping + inputs.fastq1_name + "_val_1.fq"
                self.tg_fastq2 = dirs.mapping + inputs.fastq2_name + "_val_2.fq"
            else:
                self.tg_fastq1 = dirs.mapping + inputs.fastq1_name + "_trimmed.fq"
                self.tg_fastq2 = None
            self.bt2_index = dirs.input + inputs.ref_fasta_name
            self.bt2_alignment_output = dirs.mapping + "aligned.sam"
            self.picard_bam_output = dirs.mapping + "aligned.bam"
            self.picard_sam_output = dirs.mapping + "converted.sam"
            self.picard_sorted_bam_output = dirs.mapping + "aligned_sorted.bam"
            self.picard_metrics_output = dirs.mapping + "metrics.txt"

    class _Map(object):
        def __init__(self):
            self.overwrite = False
            self.skip = False
            self.skip_fastqc = False
            self.skip_trim_galore = False
            self.tg_q_cutoff = 0
            self.bt2_alignment_args = (
                "--local;--no-unal;--no-discordant;--no-mixed;-X 1000;-L 12;-p 16"
            )
            self.description = {
                'overwrite'         : 'overwrite mapping calculation',
                'skip'              : 'do not perform sequence mapping, not recommended',
                'skip_fastqc'       : 'do not run fastqc for quality control of sequence data',
                'skip_trim_galore'  : 'do not run trim galore to remove adapter sequences at ends',
                'bt2_alignment_args': 'TODO',
                'tg_q_cutoff'       : 'TODO',
            }

    class _BitVector(object):
        def __init__(self):
            self.skip = False
            self.overwrite = False
            self.qscore_cutoff = 25
            self.num_of_surbases = 10
            self.map_score_cutoff = 15
            self.mutation_count_cutoff = 10
            self.percent_length_cutoff = 0.01
            self.summary_output_only = False
            self.plot_sequence = False
            self.miss_info = "."
            self.ambig_info = "?"
            self.nomut_bit = "0"
            self.del_bit = "1"

            self.description = {
                'overwrite'            : "overwrite bit vector calculation",
                'skip'                 : "skip bit vector generation step, not recommended",
                'qscore_cutoff'        : "quality score of read nucleotide, sets to ambigious if under this val",
                'num_of_surbases'      : "number of bases around a mutation",
                'map_score_cutoff'     : "map alignment score cutoff for a read, read is discarded if under this value",
                "mutation_count_cutoff": "maximum number of mutations in a read allowable",
                "percent_length_cutoff": "read is discarded if less than this percent of a ref sequence is included",
                "plot_sequence"        : "",
                "summary_output_only"  : ""
            }

    class _EMCluster(object):
        def __init__(self):
            pass

    def __parse_param_file(self, pf_path, p):
        if not os.path.isfile(pf_path):
            log_error_and_exit("parameter file {} supplied does not exist!".format(pf_path))
        f = open(pf_path)
        params = yaml.load(f, Loader=yaml.FullLoader)
        for group, vals in params.items():
            for param in vals:
                k, v = next(iter(param.items()))
                p.update(group + "." + k, v)

    def __parse_other_cli(self, args, p):
        # mapping
        if args['bt2_alignment_args'] is not None:
            p.map.bt2_alignment_args = args['bt2_alignment_args']
        p.map.skip_fastqc = args['skip_fastqc']
        p.map.skip_trim_galore = args['skip_trim_galore']
        # bit vector
        if args['qscore_cutoff'] is not None:
            p.bit_vector.qscore_cutoff = int(args['qscore_cutoff'])
        if args['num_of_surbases'] is not None:
            p.bit_vector.num_of_surbases = int(args['num_of_surbases'])
        if args['map_score_cutoff'] is not None:
            p.bit_vector.map_score_cutoff = int(args['map_score_cutoff'])
        if args['percent_length_cutoff'] is not None:
            p.bit_vector.percent_length_cutoff = float(args['percent_length_cutoff'])
        if args['mutation_count_cutoff'] is not None:
            p.bit_vector.mutation_count_cutoff = int(args['mutation_count_cutoff'])
        if args['plot_sequence']:
            p.bit_vector.plot_sequence = True
        if args['summary_output_only']:
            p.bit_vector.summary_output_only = True

    def get_parameters(self, args):
        input_files = validate_files(args)
        inputs = ParametersFactory._Inputs(input_files)
        if 'dot_bracket' in args:
            inputs.dot_bracket = args["dot_bracket"]
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
            p.map.overwrite = True
            p.bit_vector.overwrite = True
        if args['bv_overwrite']:
            p.bit_vector.overwrite = True
        if args['param_file'] is not None:
            self.__parse_param_file(args['param_file'], p)
        if args['log_level'] is not None:
            p.log_level = logger.str_to_log_level(args['log_level'])
        p.restore_org_behavior = args['restore_org_behavior']
        self.__parse_other_cli(args, p)
        return p


class Parameters(object):
    def __init__(
            self,
            inputs: ParametersFactory._Inputs,
            dirs: ParametersFactory._Dirs,
            files: ParametersFactory._Files,
    ):
        # general global options
        self.restore_org_behavior = False
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
        valid = set("map,files,dirs,bit_vector".split(","))
        if sub not in valid:
            log_error_and_exit(log, "unknown parameter group: {}".format(sub))
        sub_obj = self.__dict__[sub]
        if param not in sub_obj.__dict__:
            log_error_and_exit(
                    log, "unknown parameter {} in group {}".format(param, sub)
            )
        else:
            sub_obj.__dict__.update({param: val})
            log.info(
                    "setting parameter {}.{} to {} specified by parameter file".format(
                            sub, param, val
                    )
            )

    def to_yaml_file(self, fname):
        output = {}
        output["dirs"] = self.__get_args(self.dirs, "resources".split(","))
        output["map"] = self.__get_args(self.map, "description".split(","))
        output["bit_vector"] = self.__get_args(self.bit_vector,
                                               "description,miss_info,ambig_info,nomut_bit,del_bit".split(","))
        f = open(fname, "w")
        yaml.dump(output, f)
        f.close()

    def __get_args(self, group, skip):
        args = []
        for k, v in group.__dict__.items():
            if k in skip:
                continue
            args.append({k: v})
        return args


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
