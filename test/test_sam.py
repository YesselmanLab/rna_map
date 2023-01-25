"""
test sam module
"""
import os
from rna_map.sam import SingleSamIterator
from rna_map.util import fasta_to_dict

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def test_single_sam_iterator():
    fa_path = os.path.join(TEST_DIR, "resources", "case_1", "test.fasta")
    ref_seqs = fasta_to_dict(fa_path)
    sam_path = os.path.join(
        TEST_DIR,
        "resources",
        "case_1",
        "output",
        "Mapping_Files",
        "aligned.sam",
    )
    sam_iter = SingleSamIterator(sam_path, ref_seqs)
    read = next(sam_iter)[0]
    assert read.qname == 'FS10000899:22:BPG61606-0731:1:1101:1200:1000'
    assert read.flag == '0'
    assert read.rname == 'mttr-6-alt-h3'
    assert read.pos == 1
    assert read.mapq == 44
    assert read.cigar == '134M12S'
    assert read.rnext == '*'
    assert read.pnext == 0
    assert read.tlen == 0