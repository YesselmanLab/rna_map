from dreem import settings
from dreem import logger

import os
import re
import numpy as np

from Bio import SeqIO

log = logger.init_logger("bit_vector.py")


class MutationHistogram(object):
    def __init__(self, name, sequence, data_type):
        self.name = name
        self.sequence = sequence
        self.data_type = data_type
        self.num_reads = 0
        self.mut_bases = np.zeros(len(sequence))
        self.info_bases = np.zeros(len(sequence))

    @classmethod
    def from_file(cls, file_name):
        pass

    def to_file(self):
        pass

    def record_bit_vector(self, bit_vector):
        pass


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
    def __init__(self, qscore_cutoff, num_of_surbases):
        self._cigar_pattern = re.compile(r"(\d+)([A-Z]{1})")
        self._phred_qscores = self._parse_phred_qscore_file(
            settings.get_lib_path() + "/resources/phred_ascii.txt"
        )
        self._qscore_cutoff = qscore_cutoff
        self._num_of_surbases = num_of_surbases
        self._miss_info = "."
        self._ambig_info = "?"
        self._nomut_bit = "0"
        self._del_bit = "1"

    def _parse_phred_qscore_file(self, qscore_filename):
        phred_qscore = {}
        qscore_file = open(qscore_filename)
        qscore_file.readline()  # Ignore header line
        for line in qscore_file:
            line = line.strip().split()
            score, symbol = int(line[0]), line[1]
            phred_qscore[symbol] = score
        qscore_file.close()
        return phred_qscore

    def get_bit_vector(self, read: Mate, ref_seq: str):
        bitvector_mate = {}
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
                    if self._phred_qscores[q_scores[j]] > self._qscore_cutoff:
                        if read_seq[j] != ref_seq[i - 1]:
                            bitvector_mate[i] = read_seq[j]
                        else:
                            bitvector_mate[i] = self._nomut_bit
                    else:
                        bitvector_mate[i] = self._ambig_info
                    i += 1
                    j += 1
            elif desc == "D":  # Deletion
                for k in range(length - 1):
                    bitvector_mate[i] = self._ambig_info
                    i += 1
                is_ambig = self._calc_ambig_reads(
                    ref_seq, i, length, self._num_of_surbases
                )
                if is_ambig:
                    bitvector_mate[i] = self._ambig_info
                else:
                    bitvector_mate[i] = self._del_bit
                i += 1
            elif desc == "I":  # Insertion
                j += length  # Update read index
            elif desc == "S":  # soft clipping
                j += length  # Update read index
                if op_index == len(cigar_ops) - 1:  # Soft clipped at the end
                    for k in range(length):
                        bitvector_mate[i] = self._miss_info
                        i += 1
            else:
                log.warn("unknown cigar op encounters: {}".format(desc))
                return {}
            op_index += 1
        return bitvector_mate

    def _parse_cigar(self, cigar_string):
        return re.findall(self._cigar_pattern, cigar_string)

    def _calc_ambig_reads(self, ref_seq, i, length, num_surBases):
        orig_del_start = i - length + 1
        orig_sur_start = orig_del_start - num_surBases
        orig_sur_end = i + num_surBases
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


class PairedBitVectorGenerator(BitVectorGenerator):
    def __init__(self, qscore_cutoff, num_of_surbases):
        super().__init__(qscore_cutoff, num_of_surbases)
        self._bases = ["A", "C", "G", "T"]

    def get_paired_bit_vector(self, read_1, read_2, ref_seq):
        bit_vector_1 = self.get_bit_vector(read_1, ref_seq)
        bit_vector_2 = self.get_bit_vector(read_2, ref_seq)
        bit_vector = dict(bit_vector_1)
        for pos, bit in bit_vector_2.items():
            if pos not in bit_vector:  # unique to bit_vector_2
                bit_vector[pos] = bit
            elif bit != bit_vector[pos]:  # keys in both and bits not the same
                bits = set([bit_vector_1[pos], bit])
                if self._nomut_bit in bits:  # one of the bits is not mutated take that
                    bit_vector[pos] = self._nomut_bit
                # one of the bits is ambig take the other
                elif self._ambig_info in bits:
                    other_bit = list(bits - set(self._ambig_info))[0]
                    bit_vector[pos] = other_bit
                # one of the bits is missing take the other
                elif self._miss_info in bits:
                    other_bit = list(bits - set(self._miss_info))[0]
                    bit_vector[pos] = other_bit
                # both bits are mutations and different set to "?"
                elif bit_vector_1[pos] in self._bases and bit in self._bases:
                    bit_vector[pos] = self._ambig_info
                else:
                    raise ValueError(
                        "unable to merge bit_vectors with bits: {} {}".format(
                            bit_vector_1[pos], bit
                        )
                    )

        return bit_vector


def parse_fasta_file(fasta_file):
    """
    Parse a FASTA file
    Args:
        fasta_file (string): Path to FASTA file
    Returns:
        refs_seq (dict): Sequences of the ref genomes in the file
    """

    refs_seq = {}
    with open(fasta_file, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            refs_seq[record.id] = str(record.seq)
    return refs_seq


# TODO probably can just use samtools instead
def __setup_sam_file(file_locations):
    print("Converting BAM file to SAM file format")
    picard_path = settings.get_lib_path() + "/resources/picard.jar"
    convert_cmd = "java -jar {} SamFormatConverter I={} O={} >> log.txt"
    convert_cmd = convert_cmd.format(
        picard_path, file_locations.bam_file, file_locations.new_sam_file
    )
    os.system(convert_cmd)


def get_bit_vectors(file_locations, p):
    ref_seqs = parse_fasta_file(file_locations.fasta)
    mut_histos = {}
    bit_vector_writers = {}
    for ref_name, seq in ref_seqs.items():
        mut_histos[ref_name] = MutationHistogram(ref_name, seq, "DMS")
        bit_vector_writers[ref_name] = BitVectorFileWriter(
            file_locations.bitvector_files_dir, ref_name, seq, "DMS"
        )
    # __setup_sam_file(file_locations)
    ignore_lines = len(ref_seqs.keys()) + 2
    sam_file = open(file_locations.new_sam_file, "r")
    for line_index in range(ignore_lines):  # Ignore header lines
        sam_file.readline()

    if p.paired:
        parse_bit_vectors_paired(sam_file, mut_histos, bit_vector_writers, p)
    else:
        raise ValueError("not supported")


def parse_bit_vectors_paired(sam_file, mut_histos, bit_vector_writers, p):
    bvg = PairedBitVectorGenerator(p.qscore_cutoff, p.sur_bases)
    read_1_line, read_2_line = "", ""
    while True:
        read_1_line = sam_file.readline().strip()
        read_2_line = sam_file.readline().strip()
        if len(read_1_line) == 0 or len(read_2_line) == 0:
            break
        read_1 = Mate(read_1_line)
        read_2 = Mate(read_2_line)
        # check if reads are paired
        if not (
            read_1.PNEXT == read_2.POS
            and read_1.RNAME == read_2.RNAME
            and read_1.RNEXT == "="
        ):
            log.debug("mate_2 is inconsistent with mate_1 SKIPPING!")
            continue
        if not (read_1.QNAME == read_2.QNAME and read_1.MAPQ == read_2.MAPQ):
            log.debug("mate_2 is inconsistent with mate_1 SKIPPING!")
            continue
        if read_1.RNAME not in mut_histos:
            log.error("unknown ref sequence: " + read_1.RNAME)
            exit()
        bit_vector = bvg.get_paired_bit_vector(
            read_1, read_2, mut_histos[read_1.RNAME].sequence
        )
        bit_vector_writers[read_1.RNAME].write_bit_vector(read_1.QNAME, bit_vector)
