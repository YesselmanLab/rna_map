from dreem import settings
from dreem.logger import *

log = init_logger("sam.py")


class AlignedRead:
    def __init__(self, line):
        line_split = line.split()
        if len(line_split) < 11:
            raise ValueError("cannot setup Mate object from split, its too short")
        self.QNAME = line_split[0]
        self.FLAG = line_split[1]
        self.RNAME = line_split[2]
        self.POS = int(line_split[3])
        self.MAPQ = int(line_split[4])
        self.CIGAR = line_split[5]
        self.RNEXT = line_split[6]
        self.PNEXT = int(line_split[7])
        self.TLEN = line_split[8]
        self.SEQ = line_split[9]
        self.QUAL = line_split[10]
        self.MDSTRING = line_split[11].split(":")[2]


# TODO add stats to print at the end
class SamIterator(object):
    def __init__(self, samfile_path):
        self._f = open(samfile_path)
        self._good = True

    def get_next(self):
        pass

    def is_good(self):
        return self._good


class SingleSamIterator(object):
    def __init__(self, samfile_path, ref_seqs):
        self._f = open(samfile_path)
        ignore_lines = len(ref_seqs.keys()) + 2
        for line_index in range(ignore_lines):  # Ignore header lines
            self._f.readline()
        self._read_1_line = ""
        self._read_1 = None

    def __iter__(self):
        return self

    def __next__(self):
        self._read_1_line = self._f.readline().strip()
        if len(self._read_1_line) == 0:
            raise StopIteration
        self._read_1 = AlignedRead(self._read_1_line)
        return [self._read_1]

class PairedSamIterator(object):
    def __init__(self, samfile_path, ref_seqs):
        self._f = open(samfile_path)
        ignore_lines = len(ref_seqs.keys()) + 2
        for line_index in range(ignore_lines):  # Ignore header lines
            self._f.readline()
        self._read_1_line = ""
        self._read_2_line = ""
        self._read_1 = None
        self._read_2 = None

    def __iter__(self):
        return self

    def __next__(self):
        self._read_1_line = self._f.readline().strip()
        self._read_2_line = self._f.readline().strip()
        if len(self._read_1_line) == 0 or len(self._read_2_line) == 0:
            raise StopIteration
        self._read_1 = AlignedRead(self._read_1_line)
        self._read_2 = AlignedRead(self._read_2_line)
        # check if reads are paired
        if not (
            self._read_1.PNEXT == self._read_2.POS
            and self._read_1.RNAME == self._read_2.RNAME
            and self._read_1.RNEXT == "="
        ):
            log.warning(
                "mate_2 is inconsistent with mate_1 for read: {} SKIPPING!".format(
                    self._read_1.QNAME
                )
            )
            self.__next__()
        if not (
            self._read_1.QNAME == self._read_2.QNAME
            and self._read_1.MAPQ == self._read_2.MAPQ
        ):
            log.warning(
                "mate_2 is inconsistent with mate_1 for read: {} SKIPPING!".format(
                    self._read_1.QNAME
                )
            )
            self.__next__()
        return [self._read_1, self._read_2]
