import os
import re
import numpy as np
import pickle
import random
import datetime

import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
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
        self.structure = None
        self.data_type = data_type
        self.num_reads = 0
        self.mut_bases = np.zeros(len(sequence) + 1)
        self.info_bases = np.zeros(len(sequence) + 1)
        self.del_bases = np.zeros(len(sequence) + 1)
        self.ins_bases = np.zeros(len(sequence) + 1)
        self.cov_bases = np.zeros(len(sequence) + 1)
        self.mod_bases = {
            "A": np.zeros(len(sequence)),
            "C": np.zeros(len(sequence)),
            "G": np.zeros(len(sequence)),
            "T": np.zeros(len(sequence)),
        }
        self.start = start
        self.end = end
        if self.start is None:
            self.start = 1
        if self.end is None:
            self.end = len(self.sequence)

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
    xaxis_coordinates = np.array(xaxis_coordinates)
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


def plot_modified_bases(mh: MutationHistogram, p: Parameters):
    xaxis_coordinates = [i for i in range(mh.start, mh.end + 1)]
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )
    modbases_data = []
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    for base in ["A", "C", "G", "T"]:
        y_list = [mh.mod_bases[base][pos] for pos in range(mh.start, mh.end)]
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


def plot_mutation_histogram(mh: MutationHistogram, p: Parameters):
    xaxis_coordinates = [i for i in range(mh.start, mh.end + 1)]
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )
    mut_hist_data = [go.Histogram(x=mh.mut_bases)]
    mut_hist_layout = go.Layout(
        title="Mutations: " + mh.name,
        xaxis=dict(title="Number of mutations per read"),
        yaxis=dict(title="Abundance"),
    )
    mut_hist_fig = go.Figure(data=mut_hist_data, layout=mut_hist_layout)
    plotly.offline.plot(
        mut_hist_fig,
        filename=file_base_name + "mutation_histogram.html",
        auto_open=False,
    )


def plot_population_avg(mh: MutationHistogram, p: Parameters):
    xaxis_coordinates = [i for i in range(mh.start, mh.end + 1)]
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )
    popavg_filename = file_base_name + "popavg_reacts.txt"
    popavg_file = open(popavg_filename, "w")
    popavg_file.write("Position\tMismatches\tMismatches + Deletions\n")
    delmut_y, mut_y = [], []
    for pos in range(mh.start, mh.end + 1):
        try:
            delmut_frac = (mh.del_bases[pos] + mh.mut_bases[pos]) / mh.info_bases[pos]
            mut_frac = mh.mut_bases[pos] / mh.info_bases[pos]
        except ZeroDivisionError:
            delmut_frac = 0.0
            mut_frac = 0.0
        delmut_y.append(delmut_frac)
        mut_y.append(mut_frac)
        mut_frac, delmut_frac = round(mut_frac, 5), round(delmut_frac, 5)
        s = "{}\t{}\t{}\n".format(pos, mut_frac, delmut_frac)
        popavg_file.write(s)
    popavg_file.close()

    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    colors = []
    ref_bases = []
    for i in range(mh.start, mh.end + 1):
        if i >= len(mh.sequence):
            continue
        colors.append(cmap[mh.sequence[i - 1]])
        ref_bases.append(mh.sequence[i - 1])
    delmut_trace = go.Bar(
        x=xaxis_coordinates,
        y=delmut_y,
        text=ref_bases,
        marker=dict(color=colors),
        showlegend=False,
    )
    mut_trace = go.Bar(
        x=xaxis_coordinates,
        y=mut_y,
        text=ref_bases,
        marker=dict(color=colors),
        showlegend=False,
    )
    title1 = "Mismatches + Deletions: " + mh.name
    title2 = "Mismatches: " + mh.name
    mut_fig = make_subplots(rows=2, cols=1, subplot_titles=(title1, title2))
    mut_fig.append_trace(delmut_trace, 1, 1)
    mut_fig.append_trace(mut_trace, 2, 1)
    max_y = max(mut_y + delmut_y)
    mut_fig["layout"]["xaxis1"].update(title="Position")
    mut_fig["layout"]["xaxis2"].update(title="Position")
    mut_fig["layout"]["yaxis1"].update(title="Fraction", range=[0, max_y])
    mut_fig["layout"]["yaxis2"].update(title="Fraction", range=[0, max_y])
    seqs = list(mh.sequence)
    if mh.structure is not None:
        db  = list(mh.structure)
    else:
        db = " "*len(seqs)
    mut_fig.update_xaxes(
            tickvals=xaxis_coordinates,
            ticktext=["%s<br>%s" % (x, y) for (x, y) in zip(seqs, db)],
            tickangle = 0
    )
    plotly.offline.plot(
        mut_fig, filename=file_base_name + "pop_avg.html", auto_open=False,
    )


# analysis functions ###############################################################


def generate_quality_control_file(mh: MutationHistogram, p: Parameters):
    file_base_name = (
        p.dirs.bitvector + mh.name + "_" + str(mh.start) + "_" + str(mh.end) + "_"
    )
    qc_filename = file_base_name + "Quality_Control.txt"
    qc_file = open(qc_filename, "w")

    # Read coverage
    qc_file.write(mh.name + " has " + str(mh.num_reads) + " reads mapping to it")
    qc_file.write(". This is: ")
    if mh.num_reads < 50000:
        qc_file.write("BAD.\n")
    elif 50000 <= mh.num_reads < 100000:
        qc_file.write("MEDIUM.\n")
    else:
        qc_file.write("GOOD.\n")

    # Signal-noise ratio
    A_frac = mh.sequence.count("A") / len(mh.sequence)
    G_frac = mh.sequence.count("G") / len(mh.sequence)
    T_frac = mh.sequence.count("T") / len(mh.sequence)
    C_frac = mh.sequence.count("C") / len(mh.sequence)
    sig = (A_frac * np.sum(mh.mod_bases["A"])) + (C_frac * np.sum(mh.mod_bases["C"]))
    noise = (T_frac * np.sum(mh.mod_bases["T"])) + (G_frac * np.sum(mh.mod_bases["G"]))
    denom = sig + noise
    sig_noise = round(sig / denom, 2)
    qc_file.write("The signal-to-noise ratio for the sample is: " + str(sig_noise))
    qc_file.write(". This is: ")
    if sig_noise < 0.75:
        qc_file.write("BAD.\n")
    elif 0.75 <= sig_noise < 0.9:
        qc_file.write("MEDIUM.\n")
    else:
        qc_file.write("GOOD.\n")

    # Distribution of coverage
    qc_file.write("Distribution of coverage:\n")
    m = max(mh.cov_bases)
    norm_read_cov = [i / m for i in mh.cov_bases]
    n1, n2, n3 = 0, 0, 0
    for cov in norm_read_cov:
        if cov < 0.5:
            n1 += 1
        elif 0.5 <= cov < 0.75:
            n2 += 1
        elif 0.75 <= cov:
            n3 += 1
    n1 = str(round(n1 * 100 / len(mh.cov_bases), 2))
    n2 = str(round(n2 * 100 / len(mh.cov_bases), 2))
    n3 = str(round(n3 * 100 / len(mh.cov_bases), 2))
    n1_s = "{}% of bases have less than {}% of reads mapping to them\n"
    n1_s = n1_s.format(n1, 50)
    n2_s = "{}% of bases have between {}% and {}% of reads mapping to them\n"
    n2_s = n2_s.format(n2, 50, 75)
    n3_s = "{}% of bases have greater than {}% of reads mapping to them\n\n"
    n3_s = n3_s.format(n3, 75)
    qc_file.write(n1_s)
    qc_file.write(n2_s)
    qc_file.write(n3_s)

    # Info on numbers:
    qc_file.write("FOR REFERENCE:\n")

    # Read coverage
    qc_file.write("Read coverage:\n")
    qc_file.write("Number of reads < 50000: BAD\n")
    qc_file.write("50000 < Number of reads < 100000: MEDIUM\n")
    qc_file.write("Number of reads > 100000: GOOD\n\n")

    # Signal-to-noise ratio
    qc_file.write("Signal-to-noise ratio:\n")
    qc_file.write("Signal-noise ratio < 0.75: BAD\n")
    qc_file.write("0.75 < Signal-noise ratio < 0.9: MEDIUM\n")
    qc_file.write("Signal-noise ratio > 0.9: GOOD\n\n")

    qc_file.write(
        "If you are only interested in the population average"
        + " and not clustering, 1000-10000 reads might be sufficient.\n\n"
    )

    qc_file.close()


class BitVectorFileWriter(object):
    def __init__(self, path, name, sequence, data_type, start, end):
        self.start = start
        self.end = end
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
        for pos in range(self.start, self.end + 1):
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
            settings.get_py_path() + "/resources/phred_ascii.txt"
        )
        self.__bases = ["A", "C", "G", "T"]

    # TODO not big on this ... streamline somehow?
    def __setup_params(self, p: Parameters):
        self.__csv = p.ins.csv
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
        self.__generate_all_bit_vectors()
        for mh in self._mut_histos.values():
            plot_read_coverage(mh, p)
            plot_modified_bases(mh, p)
            plot_mutation_histogram(mh, p)
            plot_population_avg(mh, p)
            generate_quality_control_file(mh, p)

    def __generate_all_bit_vectors(self):
        self._mut_histos = {}
        bit_vector_pickle_file = self._p.dirs.bitvector + "mutation_histos.p"
        if os.path.isfile(bit_vector_pickle_file) and not self._p.overwrite:
            log.info(
                "SKIPPING bit vector generation, it has run already! specify -overwrite "
                + "to rerun"
            )
            with open(bit_vector_pickle_file, "rb") as handle:
                self._mut_histos = pickle.load(handle)
            return

        self._bit_vector_writers = {}
        for ref_name, seq in self._ref_seqs.items():
            self._mut_histos[ref_name] = MutationHistogram(
                ref_name, seq, "DMS", 1, len(seq)
            )
            self._bit_vector_writers[ref_name] = BitVectorFileWriter(
                self._p.dirs.bitvector, ref_name, seq, "DMS", 1, len(seq)
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

        if self.__csv is not None:
            df = pd.read_csv(self.__csv)
            for i, row in df.iterrows():
                if row["name"] in self._mut_histos:
                    self._mut_histos[row["name"]].structure = row["structure"]
        f = open(self._p.dirs.bitvector + "mutation_histos.p", "wb")
        pickle.dump(self._mut_histos, f)

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
