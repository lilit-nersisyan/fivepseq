"""
This module is meant to keep properties of annotation files and functions associated with those in one place .
"""
from preconditions import preconditions
from pysam import AlignmentFile
from plastid import BAMGenomeArray
from plastid import FivePrimeMapFactory
from fivepseq import config


class Alignment:
    alignment_file = None
    bam_array = None

    @preconditions(lambda alignment_file: isinstance(alignment_file, AlignmentFile))
    def __init__(self, alignment_file):
        """
        Initiates an Alignment class object with the given pysam.AlignmentFile.
        Creates and stores plastid BAMGenomeArray from the alignment_file, using a plastid fivePrimeMapping factory.
        :param alignment_file: pysam.AlignmentFile
        """
        self.alignment_file = alignment_file
        self.bam_array = BAMGenomeArray(alignment_file, mapping=FivePrimeMapFactory())
        print alignment_file.filename
        config.logger.debug("Read in the alignment file %s with %d chromosomes" % (alignment_file.filename, len(self.bam_array.chroms())))
