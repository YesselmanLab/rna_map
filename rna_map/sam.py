"""
custom module for parsing sam files probably should use existing package
"""

from dataclasses import dataclass
from rna_map import logger

log = logger.get_logger("SAM")


@dataclass(frozen=True, order=True)
class AlignedRead:
    qname: str
    flag: str
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str
    qual: str
    md_string: str


def get_aligned_read_from_line(line):
    """
    Get an AlignedRead object from a line of a sam file
    """
    spl = line.split()
    if len(spl) < 11:
        raise ValueError(
            "cannot setup AlignRead object from split, its too short"
        )
    return AlignedRead(
        spl[0],
        spl[1],
        spl[2],
        int(spl[3]),
        int(spl[4]),
        spl[5],
        spl[6],
        int(spl[7]),
        int(spl[8]),
        spl[9],
        spl[10],
        spl[11].split(":")[2],
    )


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
        self._read_1 = get_aligned_read_from_line(self._read_1_line)
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
        self._read_1 = get_aligned_read_from_line(self._read_1_line)
        self._read_2 = get_aligned_read_from_line(self._read_2_line)

        # check if reads are paired
        if not (
            self._read_1.pnext == self._read_2.pos
            and self._read_1.rname == self._read_2.rname
            and self._read_1.rnext == "="
        ):
            log.warning(
                f"mate_2 is inconsistent with mate_1 for read: "
                f"{self._read_1.qname} SKIPPING!"
            )
            self.__next__()
        if not (
            self._read_1.qname == self._read_2.qname
            and self._read_1.mapq == self._read_2.mapq
        ):
            log.warning(
                f"mate_2 is inconsistent with mate_1 for read: "
                f"{self._read_1.qname} SKIPPING!"
            )
            self.__next__()
        return [self._read_1, self._read_2]
