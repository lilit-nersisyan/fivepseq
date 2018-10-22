"""
The fivepseq entry point.
"""

# PORT
# argparse was written for argparse in Python 3.
# A few details are different in 2.x, especially some exception messages, which were improved in 3.x.
import argparse
import numpy
import os

import logging
import shutil
import time

import preconditions

import config
from fivepseq.util.formatting import pad_spaces
from logic.structures.counts import FivePSeqCounts
from util.readers import BamReader, AnnotationReader


class FivepseqArguments:
    """
    This class is for reading, parsing, storing and printing command line arguments.
    """

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
        parser.add_argument("--span",
                            help="the number of bases to span around a genomic position",
                            type = int,
                            required=False,
                            default=20)
        parser.add_argument("-s", "-geneset",
                            help="the file containing gene names of interest",
                            type=str,
                            required=False)
        parser.add_argument("--force",
                            help="silently overwrite files in the output directory",
                            action="store_true",
                            default=False,
                            required=False)
        parser.add_argument("--log",
                            help="set logging level",
                            choices=["INFO", "DEBUG", "info", "debug"],
                            default="INFO",
                            required=False)

        config.args = parser.parse_args()
        config.bam = os.path.abspath(config.args.b)
        config.genome = os.path.abspath(config.args.g)
        config.annot = os.path.abspath(config.args.a)
        config.span_size = config.args.span
        config.out_dir = os.path.abspath(config.args.o)

    def print_args(self):
        """
        Prints the command line arguments in a formatted manner.
        :return:
        """
        # NOTE The files are not checked for existence at this stage.
        # NOTE The absolute_path function simply merges the relative paths with the working directory.
        print "%s%s" % (pad_spaces("\tAlignment file:"), os.path.abspath(config.args.b))
        print "%s%s" % (pad_spaces("\tGenome file:"), os.path.abspath(config.args.g))
        print "%s%s" % (pad_spaces("\tAnnotation file:"), os.path.abspath(config.args.a))

        print "%s%s" % (pad_spaces("\tSpan size:"), config.span_size)
        if config.args.s is not None:
            print "%s%s" % (pad_spaces("\tGene set file:"), os.path.abspath(config.args.s))
        print "%s%s" % (pad_spaces("\tOutput directory:"), os.path.abspath(config.args.o))

        print "%s%s" % (pad_spaces("\tLogging level:"), config.args.log.upper())


def setup_output_dir():
    """
    This method attempts to set the output directory for fivepseq reports.

    If the user-specified output_dir does not exist and the parent directory does, then the output will be stored there.

    If the user-specified output_dir does not exist and the parent directory also, then an IOError is raised.

    If the user-specified output_dir exists and the force_overwrite (f) option is set to true,
    attempts to remove all the content of that directory. Exception is raised in case of failure.

    If the user specified output_dir exists and the force_overwrite (f) option is false,
    appends "+" signs to the  end of the directory name recursively, until no such directory exists.
    Creates a new directory and stores the output there.
    :param output_dir:
    :return:
    """
    output_dir = config.out_dir
    if os.path.exists(output_dir):
        if config.args.force:
            print ("\nWARNING:\tThe output directory %s already exists.\n\tFiles will be overwritten." % output_dir)
            try:
                shutil.rmtree(output_dir)
            except Exception as e:
                error_message = "ERROR:\tCould not remove directory %s. Reason:%s" % (output_dir, e.message)
                raise Exception(error_message)
        else:
            while os.path.exists(output_dir):
                output_dir += "+"
            print "\nWARNING:\tOutput directory %s already exists.\n\tFiles will be stored in %s" % (
                config.out_dir, output_dir)
            config.out_dir = output_dir

    try:
        os.makedirs(output_dir)
        print "\n\tSETUP:\tSuccess"
        print "\tSETUP:\t%s%s" % (pad_spaces("Output directory:"), output_dir)
    except Exception as e:
        error_message = "\tSETUP:\tCould not create output directory %s. Reason: %s" % (output_dir, e.message)
        raise Exception(error_message)


def setup_logger():
    """
    Sets up the logger and the handlers for different logging levels.
    Set the logger level according to the argument --log passed to the argparse. Raises a ValueError in case of an invalid argument.
    Raises an exception if the files for handlers cannot be generated.
    :return:
    """

    log_file = os.path.join(config.out_dir, "fivepseq.log")
    # config.logger = logging.getLogger()
    log_level = getattr(logging, config.args.log.upper(), None)
    if not isinstance(log_level, int):
        raise ValueError('Invalid log level: %s' % config.args.log)
    # logging.basicConfig(level=log_level, filemode='w', format = '%(levelname)s:%(asctime)s:\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S')

    # INFO level file handler
    file_handler = logging.FileHandler(log_file)
    config.logger.addHandler(file_handler)
    file_handler.setLevel(logging.INFO)
    if not os.path.exists(log_file):
        raise Exception("Could not instantiate the logger. Exiting.")
    print "\tSETUP:\t%s%s" % (pad_spaces("Log file:"), log_file)

    if log_level == logging.DEBUG:
        debug_file = os.path.join(config.out_dir, "fivepseq.debug.log")
        debug_handler = logging.FileHandler(debug_file)
        debug_handler.setLevel(logging.DEBUG)
        config.logger.addHandler(debug_handler)
        if not os.path.exists(debug_file):
            raise Exception("Could not instantiate debug logger. Exiting")
        print "\tSETUP:\t%s%s" % (pad_spaces("Debug file:"), debug_file)

    print ""


def main():
    # argument handling
    fivepseq_arguments = FivepseqArguments()
    print "Call to fivepseq with arguments:"
    fivepseq_arguments.print_args()

    # startup
    setup_output_dir()
    setup_logger()
    config.logger.info("Fivepseq started")
    start_time = time.clock()

    # body
    # TODO call the respective function
    bam_reader = BamReader(config.bam)
    annotation_reader = AnnotationReader(config.annot)
    fivepseq_counts =  FivePSeqCounts(bam_reader.alignment, annotation_reader.annotation)
    fivepseq_counts.get_counts(config.span_size,FivePSeqCounts.TERM)
    fivepseq_counts.get_counts(config.span_size,FivePSeqCounts.TERM)


    # wrap-up
    elapsed_time = time.clock() - start_time
    config.logger.info("SUCCESS! Fivepseq finished in\t%s\tseconds. The report files maybe accessed at:\t\t%s "
                       % (elapsed_time, config.out_dir))

    # TODO handle KeyboardInterrupt's (try, wait, except, pass?)


if __name__ == '__main__':
    main()
