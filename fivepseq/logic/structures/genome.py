"""
This module keeps properties of genome fasta files and functions associated with those.
"""

import logging
import gzip
from Bio import SeqIO

from fivepseq import config


COMPRESSION_GZ = ".gz"
COMPRESSION_BZ = ".bz"
COMPRESSION_None = ""

class Genome:
    fasta_file = None
    genome_dict = None
    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

    def __init__(self, fasta_file, compression = COMPRESSION_None):
        """
        Initiates a genome instance from the given path to fasta file.

        :param fasta_file: str: path to fasta file
        """

        self.fasta_file = fasta_file
        if compression == COMPRESSION_GZ:
            file_stream = gzip.open(fasta_file, 'rt', encoding='utf-8')
            self.genome_dict = SeqIO.to_dict(SeqIO.parse(file_stream, "fasta"))
        else:
            self.genome_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file), "fasta"))

        if self.genome_dict is not None:
            self.logger.debug("Genome object successfully created from file %s" % fasta_file)
        else:
            error_message = "Unknown problem occurred when reading fasta file %s" % fasta_file
            self.logger.error(error_message)
            raise Exception(error_message)
