import unittest

from fivepseq.util.readers import FastaReader
from fivepseq.tests import test_files


class GenomeTest(unittest.TestCase):
    genome = None

    def setUp(self):
        fasta_reader = FastaReader(test_files.FASTA_FILE_PATH_VALID_F)
        self.genome = fasta_reader.genome

    def test_init(self):
        self.assertIsNotNone(self.genome.genome_dict)
