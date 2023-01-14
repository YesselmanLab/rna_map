"""
test bit_vector module
"""
import os
from dreem.util import fasta_to_dict
from dreem.bit_vector import BitVectorIterator, BitVectorGenerator

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def _test_bit_vector_iterator():
    """
    test bit vector iterator
    """
    fa_path = os.path.join(TEST_DIR, "resources", "case_1", "test.fasta")
    ref_seqs = fasta_to_dict(fa_path)
    sam_path = os.path.join(
        TEST_DIR,
        "resources",
        "case_1",
        "output",
        "Mapping_Files",
        "converted.sam",
    )
    bit_vector_iter = BitVectorIterator(sam_path, ref_seqs, False)
    bit_vector = next(bit_vector_iter)
    assert len(bit_vector.data) == 146


def _test_bit_vector_generator():
    """
    test bit vector generation
    """
    fa_path = os.path.join(TEST_DIR, "resources", "case_1", "test.fasta")
    sam_path = os.path.join(
        TEST_DIR,
        "resources",
        "case_1",
        "output",
        "Mapping_Files",
        "converted.sam",
    )
    bv_gen = BitVectorGenerator()
    bv_gen.run(sam_path, fa_path, False, {})
