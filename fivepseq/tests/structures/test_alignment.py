import unittest

from fivepseq.util.readers import BamReader
from fivepseq.tests import test_files


class AlignmentTest(unittest.TestCase):
    alignment = None

    def setUp(self):
        alignment_reader = BamReader(test_files.ALIGNMENT_FILE_PATH_VALID_F)
        self.alignment = alignment_reader.alignment

    def test_init(self):
        self.assertIsNotNone(self.alignment.bam_array)

    def test_num_chromosomes(self):
        self.assertTrue(self.alignment.num_chromosomes > 0)
