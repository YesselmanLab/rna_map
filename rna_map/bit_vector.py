import os
import re
import numpy as np
import pickle
from dataclasses import dataclass
from typing import Dict, List

import pandas as pd
from tabulate import tabulate

from rna_map import settings
from rna_map.mutation_histogram import (
    MutationHistogram,
    get_dataframe,
    plot_modified_bases,
    plot_mutation_histogram,
    plot_population_avg,
    plot_read_coverage,
)
from rna_map.logger import get_logger
from rna_map.sam import AlignedRead, SingleSamIterator, PairedSamIterator
from rna_map.util import parse_phred_qscore_file, fasta_to_dict

log = get_logger("BIT_VECTOR")


@dataclass(frozen=True, order=True)
class BitVector:
    reads: List[AlignedRead]
    data: Dict


@dataclass(frozen=True, order=True)
class BitVectorSymbols:
    miss_info: str = "*"
    ambig_info: str = "?"
    nomut_bit: str = "0"
    del_bit: str = "1"


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


class BitVectorIterator(object):
    """
    A class to generate bit vectors from a SAM files. Does mininal checking
    to the bitvector is acceptable.
    """

    def __init__(
        self, sam_path, ref_seqs, paired, qscore_cutoff=25, num_of_surbases=10
    ):
        if paired:
            self.__sam_iterator = PairedSamIterator(sam_path, ref_seqs)
        else:
            self.__sam_iterator = SingleSamIterator(sam_path, ref_seqs)
        self.count = 0
        self.rejected = 0
        self.__ref_seqs = ref_seqs
        self.__paired = paired
        self.__cigar_pattern = re.compile(r"(\d+)([A-Z]{1})")
        self.__phred_qscores = parse_phred_qscore_file(
            settings.get_py_path() + "/resources/phred_ascii.txt"
        )
        # params
        self.__bases = ["A", "C", "G", "T"]
        self.__qscore_cutoff = qscore_cutoff
        self.__num_of_surbases = num_of_surbases
        self.__bts = BitVectorSymbols()

    def __iter__(self):
        return self

    def __next__(self):
        self.count += 1
        reads = next(self.__sam_iterator)
        for read in reads:
            if read.rname not in self.__ref_seqs:
                raise ValueError(
                    f"read {read.qname} aligned to {read.rname} which is not in "
                    f"the reference fasta"
                )
        if self.__paired:
            data = self.__get_bit_vector_paired(reads[0], reads[1])
        else:
            data = self.__get_bit_vector_single(reads[0])
        return BitVector(reads, data)

    def __get_bit_vector_single(self, read):
        ref_seq = self.__ref_seqs[read.rname]
        bit_vector = self.__convert_read_to_bit_vector(read, ref_seq)
        return bit_vector

    def __convert_read_to_bit_vector(self, read: AlignedRead, ref_seq: str):
        bitvector = {}
        read_seq = read.seq
        q_scores = read.qual
        i = read.pos  # Pos in the ref sequence
        j = 0  # Pos in the read sequence
        cigar_ops = self._parse_cigar(read.cigar)
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
                            bitvector[i] = self.__bts.nomut_bit
                    else:
                        bitvector[i] = self.__bts.ambig_info
                    i += 1
                    j += 1
            elif desc == "D":  # Deletion
                for k in range(length - 1):
                    bitvector[i] = self.__bts.ambig_info
                    i += 1
                is_ambig = self.__calc_ambig_reads(ref_seq, i, length)
                if is_ambig:
                    bitvector[i] = self.__bts.ambig_info
                else:
                    bitvector[i] = self.__bts.del_bit
                i += 1
            elif desc == "I":  # Insertion
                j += length  # Update read index
            elif desc == "S":  # soft clipping
                j += length  # Update read index
                if op_index == len(cigar_ops) - 1:  # Soft clipped at the end
                    for k in range(length):
                        bitvector[i] = self.__bts.miss_info
                        i += 1
            else:
                log.warn("unknown cigar op encounters: {}".format(desc))
                return {}
            op_index += 1
        return bitvector

    def __get_bit_vector_paired(self, read_1, read_2):
        ref_seq = self.__ref_seqs[read_1.rname]
        bit_vector_1 = self.__convert_read_to_bit_vector(read_1, ref_seq)
        bit_vector_2 = self.__convert_read_to_bit_vector(read_2, ref_seq)
        bit_vector = self.__merge_paired_bit_vectors(bit_vector_1, bit_vector_2)
        return bit_vector

    def _parse_cigar(self, cigar_string):
        return re.findall(self.__cigar_pattern, cigar_string)

    def __calc_ambig_reads(self, ref_seq, i, length):
        orig_del_start = i - length + 1
        orig_sur_start = orig_del_start - self.__num_of_surbases
        orig_sur_end = i + self.__num_of_surbases
        orig_sur_seq = (
            ref_seq[orig_sur_start - 1 : orig_del_start - 1]
            + ref_seq[i:orig_sur_end]
        )
        for new_del_end in range(
            i - length, i + length + 1
        ):  # Alt del end points
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

    def __merge_paired_bit_vectors(self, bit_vector_1, bit_vector_2):
        bit_vector = dict(bit_vector_1)
        for pos, bit in bit_vector_2.items():
            if pos not in bit_vector:  # unique to bit_vector_2
                bit_vector[pos] = bit
            elif bit != bit_vector[pos]:  # keys in both and bits not the same
                bits = set([bit_vector_1[pos], bit])
                if (
                    self.__bts.nomut_bit in bits
                ):  # one of the bits is not mutated take that
                    bit_vector[pos] = self.__bts.nomut_bit
                # one of the bits is ambig take the other
                elif self.__bts.ambig_info in bits:
                    other_bit = list(bits - set(self.__bts.ambig_info))[0]
                    bit_vector[pos] = other_bit
                # one of the bits is missing take the other
                elif self.__bts.miss_info in bits:
                    other_bit = list(bits - set(self.__bts.miss_info))[0]
                    bit_vector[pos] = other_bit
                # both bits are mutations and different set to "?"
                elif bit_vector_1[pos] in self.__bases and bit in self.__bases:
                    bit_vector[pos] = self.__bts.ambig_info
                # mutation on one side and insertion on the other side set to "?"
                elif (
                    bit_vector_1[pos] == self.__bts.del_bit
                    and bit in self.__bases
                    or bit_vector_1[pos] in self.__bases
                    and bit == self.__bts.del_bit
                ):
                    bit_vector[pos] = self.__bts.ambig_info
                else:
                    log.warn(
                        "unable to merge bit_vectors with bits: {} {}".format(
                            bit_vector_1[pos], bit
                        )
                    )
        return bit_vector


class BitVectorGenerator(object):
    def __init__(self):
        self.__bases = ["A", "C", "G", "T"]
        self.__bts = BitVectorSymbols()

    def setup(self, params):
        self.__params = params
        self.__out_dir = os.path.join(
            params["dirs"]["output"], "BitVector_Files/"
        )
        os.makedirs(params["dirs"]["output"], exist_ok=True)
        os.makedirs(self.__out_dir, exist_ok=True)

    def run(self, sam_path, fasta, paired, csv_file):
        log.info("starting bitvector generation")
        self.__ref_seqs = fasta_to_dict(fasta)
        self.__bit_vec_iterator = BitVectorIterator(
            sam_path, self.__ref_seqs, paired
        )
        self.__mut_histos = {}
        self.__map_score_cutoff = self.__params["bit_vector"][
            "map_score_cutoff"
        ]
        self.__csv_file = csv_file
        self.__summary_only = self.__params["bit_vector"]["summary_output_only"]
        self.__rejected_out = open(self.__out_dir + "rejected_bvs.csv", "w")
        self.__rejected_out.write("qname,rname,reason,read1,read2,bitvector\n")
        # setup parameters about generating bit vectors
        self.__generate_all_bit_vectors()
        self.__generate_plots()
        self.__get_skip_summary()
        self.__write_summary_csv()

    def __write_summary_csv(self):
        cols = [
            "name",
            "reads",
            "aligned",
            "no_mut",
            "1_mut",
            "2_mut",
            "3_mut",
            "3plus_mut",
            "sn",
        ]
        df = get_dataframe(self.__mut_histos, cols)
        log.info(
            "MUTATION SUMMARY:\n"
            + tabulate(df, df.columns, tablefmt="github", showindex=False)
        )

        sum_path = os.path.join(self.__out_dir, "summary.csv")
        df.to_csv(sum_path, index=False)

    def __get_skip_summary(self):
        data = []
        cols = ["low_mapq"]
        if self.__params["stricter_bv_constraints"]:
            cols += ["short_read", "too_many_muts", "muts_too_close"]
        for mut_histo in self.__mut_histos.values():
            row = [mut_histo.name]
            for col in cols:
                row.append(mut_histo.skips[col] / mut_histo.num_reads * 100)
            data.append(row)
        df = pd.DataFrame(data, columns=["name"] + cols)
        log.info(
            "REMOVED READS:\n"
            + tabulate(df, df.columns, tablefmt="github", showindex=False)
            + "\n"
        )

    def __generate_plots(self):
        """
        Generate plots for each mutation histogram
        """
        for _, mh in self.__mut_histos.items():
            fname = f"{self.__out_dir}/{mh.name}_{mh.start}_{mh.end}_"
            if not self.__summary_only:
                df = mh.get_pop_avg_dataframe()
                plot_population_avg(df, mh.name, f"{fname}pop_avg.html")
            # TODO add generate other plots arg?
            if self.__params["restore_org_behavior"]:
                plot_modified_bases(
                    mh.get_nuc_coords(), mh.mod_bases, f"{fname}mutations.html"
                )
                plot_mutation_histogram(
                    mh.get_nuc_coords(),
                    mh.num_of_mutations,
                    f"{fname}mutation_histogram.html",
                )
                plot_read_coverage(
                    mh.get_nuc_coords(),
                    mh.get_read_coverage(),
                    f"{fname}read_coverage.html",
                )

    def __generate_all_bit_vectors(self):
        # TODO turn this overwrite check back on
        """if (
            os.path.isfile(bit_vector_pickle_file)
            and not self._p.bit_vector.overwrite
        ):
            log.info(
                "SKIPPING bit vector generation, it has run already! specify -overwrite "
                + "to rerun"
            )
            with open(bit_vector_pickle_file, "rb") as handle:
                self._mut_histos = pickle.load(handle)
            return"""

        self._bit_vector_writers = {}
        for ref_name, seq in self.__ref_seqs.items():
            self.__mut_histos[ref_name] = MutationHistogram(
                ref_name, seq, "DMS", 1, len(seq)
            )
            if not self.__summary_only:
                self._bit_vector_writers[ref_name] = BitVectorFileWriter(
                    self.__out_dir, ref_name, seq, "DMS", 1, len(seq)
                )
        if self.__csv_file != "":
            df = pd.read_csv(self.__csv_file)
            for i, row in df.iterrows():
                if row["name"] in self.__mut_histos:
                    self.__mut_histos[row["name"]].structure = row["structure"]
        for bit_vector in self.__bit_vec_iterator:
            self.__record_bit_vector(bit_vector)
        # pickle mutational histograms
        pickle_file = os.path.join(self.__out_dir, "mutation_histos.p")
        with open(pickle_file, "wb") as handle:
            pickle.dump(self.__mut_histos, handle)

    def __record_bit_vector(self, bit_vector):
        mh = self.__mut_histos[bit_vector.reads[0].rname]
        # if the reads do not meet the minimum mapping score, skip
        for read in bit_vector.reads:
            if read.mapq < self.__map_score_cutoff:
                self.__write_rejected_bit_vector(mh, bit_vector, "low_mapq")
                mh.record_skip("low_mapq")
                return
        # experimental features
        # must turn --stricter-bv-constraints to use these
        if self.__are_reads_to_short(mh, bit_vector):
            return
        if self.__too_many_mutations(mh, bit_vector):
            return
        if self.__muts_too_close(mh, bit_vector):
            return
        self.__update_mut_histo(mh, bit_vector.data)
        if not self.__params["bit_vector"]["summary_output_only"]:
            self._bit_vector_writers[bit_vector.reads[0].rname].write_bit_vector(
                bit_vector.reads[0].qname, bit_vector.data
            )

    def __update_mut_histo(self, mh: MutationHistogram, bit_vector):
        mh.num_reads += 1
        mh.num_aligned += 1
        total_muts = 0
        for pos in mh.get_nuc_coords():
            if pos not in bit_vector:
                continue
            read_bit = bit_vector[pos]
            if read_bit != self.__bts.ambig_info:
                mh.cov_bases[pos] += 1
            if read_bit in self.__bases:
                total_muts += 1
                mh.mod_bases[read_bit][pos] += 1
                mh.mut_bases[pos] += 1
            elif read_bit == self.__bts.del_bit:
                mh.del_bases[pos] += 1
            mh.info_bases[pos] += 1
        mh.num_of_mutations[total_muts] += 1

    # bit vector constraints ###################################################

    def __are_reads_to_short(self, mh, bit_vector) -> bool:
        if not self.__params["stricter_bv_constraints"]:
            return False
        cutoff = self.__params["bit_vector"]["stricter_constraints"][
            "percent_length_cutoff"
        ]
        ref_seq = self.__ref_seqs[bit_vector.reads[0].rname]
        for read in bit_vector.reads:
            per = len(read.seq) / len(ref_seq)
            if per < cutoff:
                self.__write_rejected_bit_vector(mh, bit_vector, "short_read")
                mh.record_skip("short_read")
                return True
        return False

    def __too_many_mutations(self, mh, bit_vector):
        if not self.__params["stricter_bv_constraints"]:
            return False
        cutoff = self.__params["bit_vector"]["stricter_constraints"][
            "mutation_count_cutoff"
        ]
        muts = 0
        for pos in mh.get_nuc_coords():
            if pos not in bit_vector.data:
                continue
            read_bit = bit_vector.data[pos]
            if read_bit in self.__bases:
                muts += 1
        if muts > cutoff:
            self.__write_rejected_bit_vector(mh, bit_vector, "too_many_muts")
            mh.record_skip("too_many_muts")
            return True
        return False

    def __muts_too_close(self, mh, bit_vector):
        if not self.__params["stricter_bv_constraints"]:
            return False
        cutoff = self.__params["bit_vector"]["stricter_constraints"][
            "min_mut_distance"
        ]
        for pos in range(mh.start, mh.end + 1):
            if pos not in bit_vector.data:
                continue
            read_bit = bit_vector.data[pos]
            if read_bit in self.__bases:
                for pos2 in range(pos - cutoff, pos + cutoff):
                    if pos2 == pos:
                        continue
                    if pos2 not in bit_vector.data:
                        continue
                    if bit_vector.data[pos2] in self.__bases:
                        self.__write_rejected_bit_vector(
                            mh, bit_vector, "muts_too_close"
                        )
                        mh.record_skip("muts_too_close")
                        return True
        return False

    def __write_rejected_bit_vector(self, mh, bit_vector, reason):
        read1 = bit_vector.reads[0]
        if len(bit_vector.reads) == 2:
            read2_seq = bit_vector.reads[1].seq
        else:
            read2_seq = ""
        bv_vec = []
        for nuc in mh.get_nuc_coords():
            if nuc in bit_vector.data:
                bv_vec.append(bit_vector.data[nuc])
            else:
                bv_vec.append(".")
        self.__rejected_out.write(
            f"{read1.qname},{read1.rname},{reason},{read1.seq},{read2_seq},"
            f"{''.join(bv_vec)}\n"
        )
