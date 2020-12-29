from dreem import settings


class Read:
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


class FastqParser(object):
    def __init__(self, samfile_path):
        self._f = open(samfile_path)
        self._good = True

    def get_next(self):
        pass

    def is_good(self):
        return self._good


class PairedFastqParser(FastqParser):
    pass
