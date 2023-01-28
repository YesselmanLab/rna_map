"""
test bit_vector module
"""
import os
import shutil
import pytest
from rna_map.util import fasta_to_dict
from rna_map.parameters import get_default_params
from rna_map.bit_vector import BitVectorIterator, BitVectorGenerator

TEST_DIR = os.path.dirname(os.path.realpath(__file__))



@pytest.mark.quick
def test_bit_vector_iterator():
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
        "aligned.sam",
    )
    bit_vector_iter = BitVectorIterator(sam_path, ref_seqs, False)
    bit_vector = next(bit_vector_iter)
    assert len(bit_vector.data) == 146
    count = 0
    for _ in bit_vector_iter:
        count += 1
    assert count == 2356

@pytest.mark.quick
def test_bit_vector_generator():
    """
    test bit vector generation
    """
    # setup_applevel_logger()
    fa_path = os.path.join(TEST_DIR, "resources", "case_1", "test.fasta")
    sam_path = os.path.join(
        TEST_DIR,
        "resources",
        "case_1",
        "output",
        "Mapping_Files",
        "aligned.sam",
    )
    bv_gen = BitVectorGenerator()
    bv_gen.setup(get_default_params())
    bv_gen.run(sam_path, fa_path, False, "")
    assert os.path.exists("output/BitVector_Files/summary.csv")
    shutil.rmtree("output")
