import os
import re
import numpy as np
import pickle

import plotly
import plotly.graph_objs as go
from plotly import tools
import datetime
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from dreem import settings, fastq
from dreem.parameters import Parameters
from dreem.logger import *
from dreem.util import *
from dreem.fastq import *

log = init_logger("bit_vector.py")


class MutationHistogram(object):
    def __init__(self, name, sequence, data_type, start=None, end=None):
        self.__bases = ["A", "C", "G", "T"]
        self.name = name
        self.sequence = sequence
        self.data_type = data_type
        self.num_reads = 0
        self.mut_bases = np.zeros(len(sequence))
        self.info_bases = np.zeros(len(sequence))
        self.del_bases = np.zeros(len(sequence))
        self.ins_bases = np.zeros(len(sequence))
        self.cov_bases = np.zeros(len(sequence))
        self.mod_bases = {
            "A": np.zeros(len(sequence)),
            "C": np.zeros(len(sequence)),
            "G": np.zeros(len(sequence)),
            "T": np.zeros(len(sequence)),
        }
        self.start = start
        self.end = end
        if self.start is None:
            self.start = 0
        if self.end is None:
            self.end = len(self.sequence) - 1

    @classmethod
    def from_file(cls, file_name):
        pass

    def to_file(self):
        pass

    def record_bit_vector(self, bit_vector, p):
        self.num_reads += 1
        for pos in range(self.start, self.end + 1):
            if pos not in bit_vector:
                continue
            read_bit = bit_vector[pos]
            if read_bit != p.bit_vector.ambig_info:
                self.cov_bases[pos] += 1
                self.info_bases[pos] += 1
            if read_bit in self.__bases:
                self.mod_bases[read_bit][pos] += 1
                self.mut_bases[pos] += 1
            elif read_bit == p.bit_vector.del_bit:
                self.del_bases[pos] += 1


# plotting functions ###############################################################
def plot_read_coverage(mh: MutationHistogram, p: Parameters):
    xaxis_coordinates = [i for i in range(mh.start, mh.end + 1)]
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )
    read_cov = []
    for pos in range(mh.start, mh.end + 1):
        try:
            cov_frac = mh.cov_bases[pos] / mh.num_reads
        except ZeroDivisionError:
            cov_frac = 0.0
        read_cov.append(cov_frac)
    cov_trace = go.Bar(x=xaxis_coordinates, y=read_cov)
    cov_data = [cov_trace]
    cov_layout = go.Layout(
        title="Read coverage: "
        + mh.name
        + ", Number of bit vectors: "
        + str(mh.num_reads),
        xaxis=dict(title="Position"),
        yaxis=dict(title="Coverage fraction"),
    )
    cov_fig = go.Figure(data=cov_data, layout=cov_layout)
    plotly.offline.plot(
        cov_fig, filename=file_base_name + "read_coverage.html", auto_open=False,
    )
    plotly.io.write_image(
        cov_fig, file_base_name + "read_coverage.pdf",
    )


def plot_modified_bases(mh: MutationHistogram, p: Parameters):
    xaxis_coordinates = [i for i in range(mh.start, mh.end + 1)]
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )
    modbases_data = []
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    for base in ["A", "C", "G", "T"]:
        y_list = [mh.mod_bases[base][pos] for pos in range(mh.start, mh.end + 1)]
        trace = go.Bar(
            x=xaxis_coordinates, y=y_list, name=base, marker_color=cmap[base]
        )
        modbases_data.append(trace)
    modbases_layout = go.Layout(
        title="DMS modifications: " + mh.name,
        xaxis=dict(title="Position"),
        yaxis=dict(title="Abundance"),
        barmode="stack",
    )
    modbases_fig = go.Figure(data=modbases_data, layout=modbases_layout)
    modbases_fig
    plotly.offline.plot(
        modbases_fig, filename=file_base_name + "mutations.html", auto_open=False,
    )
    # add pdf
    plotly.io.write_image(
        modbases_fig, file_base_name + "mutations.pdf",
    )


def plot_mutation_histogram():
    xaxis_coordinates = [i for i in range(mh.start, mh.end + 1)]
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )


class BitVectorFileWriter(object):
    def __init__(self, path, name, sequence, data_type):
        self.sequence = sequence
        self.f = open(path + name + "_bitvectors.txt", "w")
        self.f.write("@ref\t{}\t{}\t{}\n".format(name, sequence, data_type))
        self.f.write(
            "@coordinates:\t{},{}:{}\n".format(0, len(sequence), len(sequence))
        )
        self.f.write("Query_name\tBit_vector\tN_Mutations\n")

    def write_bit_vector(self, q_name, bit_vector):
        n_mutations = 0
        bit_string = ""
        for pos in range(0, len(self.sequence)):
            if pos not in bit_vector:
                bit_string += "."
            else:
                read_bit = bit_vector[pos]
                if read_bit.isalpha():
                    n_mutations += 1
                bit_string += read_bit
        self.f.write("{}\t{}\t{}\n".format(q_name, bit_string, n_mutations))


class BitVectorFileReader(object):
    def __init__(self):
        pass


class BitVectorGenerator(object):
    def __init__(self):
        self.__cigar_pattern = re.compile(r"(\d+)([A-Z]{1})")
        self.__phred_qscores = self.__parse_phred_qscore_file(
            settings.get_lib_path() + "/resources/phred_ascii.txt"
        )
        self.__bases = ["A", "C", "G", "T"]

    # TODO not big on this ... streamline somehow?
    def __setup_params(self, p: Parameters):
        self.__qscore_cutoff = p.bit_vector.qscore_cutoff
        self.__num_of_surbases = p.bit_vector.num_of_surbases
        self.__miss_info = p.bit_vector.miss_info
        self.__ambig_info = p.bit_vector.ambig_info
        self.__nomut_bit = p.bit_vector.nomut_bit
        self.__del_bit = p.bit_vector.del_bit
        self.__map_score_cutoff = p.bit_vector.map_score_cutoff
        self.__mutation_count_cutoff = p.bit_vector.mutation_count_cutoff
        self.__percent_length_cutoff = p.bit_vector.percent_length_cutoff

    def run(self, p: Parameters):
        log.info("starting bitvector generation")
        self._p = p
        # setup parameters about generating bit vectors
        self.__setup_params(p)
        self._ref_seqs = fasta_to_dict(self._p.ins.ref_fasta)
        log.setLevel(p.log_level)
        self.__run_picard_sam_convert()

        self._mut_histos = {}
        self._bit_vector_writers = {}
        for ref_name, seq in self._ref_seqs.items():
            self._mut_histos[ref_name] = MutationHistogram(ref_name, seq, "DMS")
            self._bit_vector_writers[ref_name] = BitVectorFileWriter(
                self._p.dirs.bitvector, ref_name, seq, "DMS"
            )
        fastq_iterator = None
        if self._p.paired:
            fastq_iterator = fastq.PairedFastqIterator(
                self._p.files.picard_sam_output, self._ref_seqs
            )
        for read in fastq_iterator:
            if self._p.paired:
                bit_vector = self.__get_bit_vector_paired(read[0], read[1])
            else:
                bit_vector = self.__get_bit_vector_single(read)
        f = open(self._p.dirs.bitvector + "mutation_histos.p", "wb")
        pickle.dump(self._mut_histos, f)

        for mh in self._mut_histos.values():
            plot_read_coverage(mh, p)
            plot_modified_bases(mh, p)

    def __get_bit_vector_single(self, read):
        raise NotImplemented()

    def __get_bit_vector_paired(self, read_1, read_2):
        if read_1.RNAME not in self._ref_seqs:
            log_error_and_exit(
                "read {} aligned to {} which is not in the reference fasta".format(
                    read_1.QNAME, read_1.RNAME
                )
            )
        ref_seq = self._ref_seqs[read_1.RNAME]
        bit_vector_1 = self.__convert_read_to_bit_vector(read_1, ref_seq)
        bit_vector_2 = self.__convert_read_to_bit_vector(read_2, ref_seq)
        bit_vector = self.__merge_paired_bit_vectors(bit_vector_1, bit_vector_2)
        self._bit_vector_writers[read_1.RNAME].write_bit_vector(
            read_1.QNAME, bit_vector
        )
        self._mut_histos[read_1.RNAME].record_bit_vector(bit_vector, self._p)
        return bit_vector

    def __run_picard_sam_convert(self):
        if os.path.isfile(self._p.files.picard_sam_output) and not self._p.overwrite:
            log.info(
                "SKIPPING picard SAM convert, it has been run already! specify "
                + "-overwrite to rerun"
            )
            return

        picard_path = self._p.dirs.resources + "/picard.jar"
        picard_cmd = (
            "java -jar "
            + picard_path
            + " SamFormatConverter I={} O={}".format(
                self._p.files.picard_bam_output, self._p.files.picard_sam_output
            )
        )
        self.__run_command("picard SAM convert", picard_cmd)

    # TODO copied from mapping ... centralize or remove
    def __run_command(self, method_name, cmd):
        log.info("running {}".format(method_name))
        log.debug(cmd)
        output, error_msg = run_command(cmd)
        if error_msg is not None:
            self.__log_error_msg_and_exit(log, method_name, error_msg)
        else:
            log.info("{} ran without errors".format(method_name))
        return output, error_msg

    def __log_error_msg_and_exit(self, log, pname, error_msg):
        log.error("{} returned error:".format(pname))
        log.error(error_msg)
        log.error("EXITING")
        exit()

    def __parse_phred_qscore_file(self, qscore_filename):
        phred_qscore = {}
        qscore_file = open(qscore_filename)
        qscore_file.readline()  # Ignore header line
        for line in qscore_file:
            line = line.strip().split()
            score, symbol = int(line[0]), line[1]
            phred_qscore[symbol] = score
        qscore_file.close()
        return phred_qscore

    def _parse_cigar(self, cigar_string):
        return re.findall(self.__cigar_pattern, cigar_string)

    def __convert_read_to_bit_vector(self, read: AlignedRead, ref_seq: str):
        bitvector = {}
        read_seq = read.SEQ
        q_scores = read.QUAL
        i = read.POS  # Pos in the ref sequence
        j = 0  # Pos in the read sequence
        cigar_ops = self._parse_cigar(read.CIGAR)
        op_index = 0
        while op_index < len(cigar_ops):
            op = cigar_ops[op_index]
            desc, length = op[1], int(op[0])
            if desc == "M":  # Match or mismatch
                for k in range(length):
                    if self.__phred_qscores[q_scores[j]] > self.__qscore_cutoff:
                        if read_seq[j] != ref_seq[i - 1]:
                            bitvector[i] = read_seq[j]
                        else:
                            bitvector[i] = self.__nomut_bit
                    else:
                        bitvector[i] = self.__ambig_info
                    i += 1
                    j += 1
            elif desc == "D":  # Deletion
                for k in range(length - 1):
                    bitvector[i] = self.__ambig_info
                    i += 1
                is_ambig = self.__calc_ambig_reads(ref_seq, i, length)
                if is_ambig:
                    bitvector[i] = self.__ambig_info
                else:
                    bitvector[i] = self.__del_bit
                i += 1
            elif desc == "I":  # Insertion
                j += length  # Update read index
            elif desc == "S":  # soft clipping
                j += length  # Update read index
                if op_index == len(cigar_ops) - 1:  # Soft clipped at the end
                    for k in range(length):
                        bitvector[i] = self.__miss_info
                        i += 1
            else:
                log.warn("unknown cigar op encounters: {}".format(desc))
                return {}
            op_index += 1
        return bitvector

    def __merge_paired_bit_vectors(self, bit_vector_1, bit_vector_2):
        bit_vector = dict(bit_vector_1)
        for pos, bit in bit_vector_2.items():
            if pos not in bit_vector:  # unique to bit_vector_2
                bit_vector[pos] = bit
            elif bit != bit_vector[pos]:  # keys in both and bits not the same
                bits = set([bit_vector_1[pos], bit])
                if self.__nomut_bit in bits:  # one of the bits is not mutated take that
                    bit_vector[pos] = self.__nomut_bit
                # one of the bits is ambig take the other
                elif self.__ambig_info in bits:
                    other_bit = list(bits - set(self.__ambig_info))[0]
                    bit_vector[pos] = other_bit
                # one of the bits is missing take the other
                elif self.__miss_info in bits:
                    other_bit = list(bits - set(self.__miss_info))[0]
                    bit_vector[pos] = other_bit
                # both bits are mutations and different set to "?"
                elif bit_vector_1[pos] in self.__bases and bit in self.__bases:
                    bit_vector[pos] = self.__ambig_info
                else:
                    raise ValueError(
                        "unable to merge bit_vectors with bits: {} {}".format(
                            bit_vector_1[pos], bit
                        )
                    )
        return bit_vector

    def __calc_ambig_reads(self, ref_seq, i, length):
        orig_del_start = i - length + 1
        orig_sur_start = orig_del_start - self.__num_of_surbases
        orig_sur_end = i + self.__num_of_surbases
        orig_sur_seq = (
            ref_seq[orig_sur_start - 1 : orig_del_start - 1] + ref_seq[i:orig_sur_end]
        )
        for new_del_end in range(i - length, i + length + 1):  # Alt del end points
            if new_del_end == i:  # Orig end point
                continue
            new_del_start = new_del_end - length + 1
            sur_seq = (
                ref_seq[orig_sur_start - 1 : new_del_start - 1]
                + ref_seq[new_del_end:orig_sur_end]
            )
            if sur_seq == orig_sur_seq:
                return True
        return False
