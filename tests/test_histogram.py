import pickle
import unittest

from dreem.bit_vector import MutationHistogram

class MutationHistogramUnittest(unittest.TestCase):

    def test(self):
        mh = MutationHistogram("test_construct", "GGGGAAAACCCC", "DMS")


def main():
    unittest.main()


if __name__ == '__main__':
    main()