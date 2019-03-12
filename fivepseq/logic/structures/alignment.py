"""
This module is meant to keep properties of alignment files and functions associated with those in one place .
"""
import logging

from plastid import BAMGenomeArray
from plastid import FivePrimeMapFactory
from preconditions import preconditions
from pysam import AlignmentFile

from fivepseq import config


class Alignment:
    alignment_file = None
    bam_array = None
    num_chromosomes = None
    logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)

    @preconditions(lambda alignment_file: isinstance(alignment_file, AlignmentFile))
    def __init__(self, alignment_file):

        """
        Initiates an Alignment class object with the given pysam.AlignmentFile.
        Creates and stores plastid BAMGenomeArray from the alignment_file, using a plastid fivePrimeMapping factory.

        :param alignment_file: pysam.AlignmentFile

        """

        self.alignment_file = alignment_file
        self.bam_array = BAMGenomeArray(alignment_file, mapping=FivePrimeMapFactory())

        # check if alignment file read successfully
        if self.bam_array is None:
            error_message = "Unknown problem occurred while reading the alignment file" % alignment_file
            self.logger.error(error_message)
            raise Exception(error_message)

        # set the number of chromosomes from the bam array
        self.num_chromosomes = len(self.bam_array.chroms())

        # report success
        self.logger.debug("Read in the alignment file %s with %d chromosomes"
                            % (alignment_file.filename, self.num_chromosomes))
