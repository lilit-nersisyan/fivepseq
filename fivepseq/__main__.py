"""
The fivepseq entry point.
"""

# PORT
# argparse was written for argparse in Python 3.
# A few details are different in 2.x, especially some exception messages, which were improved in 3.x.
import argparse
import os

import logging

from util.formatting import pad_spaces
from util.readers import BamReader


class FivepseqArguments:
    """
    This class is for reading, parsing, storing and printing command line arguments.
    """
    args = None
    bam = None
    genome = None
    annot = None
    logger = None

    def __init__(self):
        """
        Initializes FivepseqArguments with command line arguments.
        Stores those in the <args> attribute.
        """
        parser = argparse.ArgumentParser(
            description="Reports 5'-footprint properties from 5Pseq reads in a given alignment file")

        parser.add_argument("-b", "-bam",
                            help="the full path to the bam/sam file (may be compressed)",
                            type=str,
                            required=True)
        parser.add_argument("-g", "-genome",
                            help="the full path to the fa/fasta file (may be compressed)",
                            type=str,
                            required=True)
        parser.add_argument("-a", "-annotation",
                            help="the full path to the gtf/gff/gff3 file (may be compressed)",
                            type=str,
                            required=True)
        parser.add_argument("-o", "-outdir",
                            help="the output directory",
                            type=str,
                            required=False,
                            default="fivepseq_out")
        parser.add_argument("-s", "-geneset",
                            help="the file containing gene names of interest",
                            type=str,
                            required=False)
        parser.add_argument("-f", "--force_overwrite",
                            help="silently overwrite files in the output directory",
                            action="store_false",
                            required=False)
        parser.add_argument("-v", "--verbosity",
                            help="increase output verbosity",
                            action="store_false",
                            default=0)
        parser.add_argument("-lv", "--log_verbosity",
                            help="increase log verbosity",
                            type=str,
                            choices=["hi", "me", "lo"])

        self.args = parser.parse_args()
        self.bam = os.path.abspath(self.args.b)
        self.genome = os.path.abspath(self.args.g)
        self.annot = os.path.abspath(self.args.a)

    def print_args(self):
        """
        Prints the command line arguments in a formatted manner.
        :return:
        """
        # NOTE The files are not checked for existence at this stage.
        # NOTE The absolute_path function simply merges the relative paths with the working directory.
        print "%s%s" % (pad_spaces("\tAlignment file:"), os.path.abspath(self.args.b))
        print "%s%s" % (pad_spaces("\tGenome file:"), os.path.abspath(self.args.g))
        print "%s%s" % (pad_spaces("\tAnnotation file:"), os.path.abspath(self.args.a))

        if self.args.s is not None:
            print "%s%s" % (pad_spaces("\tGene set file:"), os.path.abspath(self.args.s))
        print "%s%s" % (pad_spaces("\tOutput directory:"), os.path.abspath(self.args.o))

        if self.args.verbosity:
            print "\n\tverbosity turned on"

    def logger_setup(self):
        self.logger = logging.getLogger(__name__)
        if self.args.verbosity:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.DEBUG)
        logfile = os.path.join(self.args.o, "fivepseq.log")
        handler = logging.FileHandler(logfile)
        self.logger.addHandler(handler)

def main():
    fivepseq_arguments = FivepseqArguments()
    print "Started fivepseq with arguments:"
    fivepseq_arguments.print_args()

    # TODO call the respective function
    bamreader = BamReader(fivepseq_arguments.args.b)
    print bamreader
    print "fivepseq finished"

    # TODO handle KeyboardInterrupt's (try, wait, except, pass?)


if __name__ == '__main__':
    main()
