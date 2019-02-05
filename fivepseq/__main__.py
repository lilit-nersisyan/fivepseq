"""
The fivepseq entry point.
"""
# PORT
# argparse was written for argparse in Python 3.
# A few details are different in 2.x, especially some exception messages, which were improved in 3.x.
import argparse

import os
import logging
import shutil
import time

from fivepseq import config
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts
from fivepseq.util.readers import BamReader, AnnotationReader, FastaReader
from fivepseq.util.formatting import pad_spaces
from fivepseq.logic.structures.fivepseq_counts import CountManager
from fivepseq.util.writers import FivePSeqOut


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

        # TODO parse arguments for fivepseq count and viz programs separately;
        # TODO implement main accordingly
        # parser = argparse.ArgumentParser(prog = "fivepseq")
        # sp = parser.add_subparsers()
        # sp_count = sp.add_parser('count', help = 'Computes and stores %(prog) counts')
        # sp_viz = sp.add_parser('viz', help = 'Visualizes %(prog) outputs')
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
                            type=int,
                            required=False,
                            default=100)
        parser.add_argument("-s", "-geneset",
                            help="the file containing gene names of interest",
                            type=str,
                            required=False)
        parser.add_argument("--force",
                            help="silently overwrite files in the output directory",
                            action="store_true",
                            default=False,
                            required=False)
        # parser.add_argument("--addFiles",
        #                    help="only add files if not already generated",
        #                    action="store_true",
        #                    default=False,
        #                    required=False)
        # arg_group_output = parser.add_mutually_exclusive_group()
        # arg_group_output.add_argument("--force", action="store_true")
        # arg_group_output.add_argument("--addFiles", action="store_true")

        parser.add_argument("--log",
                            help="set logging level",
                            choices=["INFO", "DEBUG", "info", "debug"],
                            default="INFO",
                            required=False)
        parser.add_argument("--ignore-cache",
                            help="read in transcript assembly and rewrite existing pickle paths",
                            action="store_true",
                            default=False,
                            required=False)

        parser.add_argument("--loci-file",
                            help="file with loci to count mapping positions relative to",
                            required=False,
                            type=str)

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
        # TODO implement this option in writers
        # elif config.args.addFiles:
        #    print ("\nWARNING:\tThe output directory %s already exists.\n\tFiles will be added if not already there." % output_dir)
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

    # FIXME probably the cache directory should reside somewhere in the fivepseq package folder
    # TODO add checks to make sure the parent directory exists and has write permissions
    # TODO (maybe should not be a problem when it is within the fivepseq package, or not?)
    config.cache_dir = os.path.join(os.path.dirname(config.out_dir), "fivepseq_cache")
    if not os.path.exists(config.cache_dir):
        os.makedirs(config.cache_dir)


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


def generate_and_store_fivepseq_counts():
    # read files
    bam_reader = BamReader(config.bam)
    annotation_reader = AnnotationReader(config.annot)  # set the break for human
    fasta_reader = FastaReader(config.genome)

    # combine objects into FivePSeqCounts object
    fivepseq_counts = FivePSeqCounts(bam_reader.alignment, annotation_reader.annotation, fasta_reader.genome,
                                     config.span_size)

    # compute and store the counts in files
    fivepseq_out = FivePSeqOut(config.out_dir)


    #   terminal counts
    term_counts = fivepseq_counts.get_count_vector_list(FivePSeqCounts.TERM)
    fivepseq_out.write_vector_list(term_counts, fivepseq_out.COUNT_TERM_FILE)
    #   start counts
    start_counts = fivepseq_counts.get_count_vector_list(FivePSeqCounts.START)
    fivepseq_out.write_vector_list(start_counts, fivepseq_out.COUNT_START_FILE)
    #   full length counts
    full_length_counts = fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH)
    fivepseq_out.write_vector_list(full_length_counts, fivepseq_out.COUNT_FULL_FILE)
    fivepseq_out.write_dict(fivepseq_counts.start_codon_dict, fivepseq_out.START_CODON_DICT_FILE)
    fivepseq_out.write_dict(fivepseq_counts.stop_codon_dict, fivepseq_out.TERM_CODON_DICT_FILE)
    fivepseq_out.write_vector(fivepseq_counts.canonical_transcript_index, fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)
    fivepseq_out.write_df_to_file(fivepseq_counts.transcript_descriptors, fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)
    #   meta counts
    fivepseq_out.write_series_to_file(fivepseq_counts.get_meta_count_series(FivePSeqCounts.TERM),
                                      fivepseq_out.META_COUNT_TERM_FILE)
    fivepseq_out.write_series_to_file(fivepseq_counts.get_meta_count_series(FivePSeqCounts.START),
                                      fivepseq_out.META_COUNT_START_FILE)
    #   frame counts
    fivepseq_out.write_df_to_file(
        CountManager.extract_count_sums_per_frame_per_transcript(full_length_counts, config.span_size,
                                                                 FivePSeqCounts.START),
        fivepseq_out.FRAME_COUNTS_START_FILE)
    fivepseq_out.write_df_to_file(
        CountManager.extract_count_sums_per_frame_per_transcript(full_length_counts, config.span_size,
                                                                 FivePSeqCounts.TERM),
        fivepseq_out.FRAME_COUNTS_TERM_FILE)


    fivepseq_out.write_df_to_file(fivepseq_counts.get_amino_acid_pauses(config.span_size, 20),
                                  fivepseq_out.AMINO_ACID_PAUSES_FILE)

    #   transcript assembly
    fivepseq_out.write_transcript_assembly_to_file(annotation_reader.annotation.transcript_assembly,
                                                   fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE)

    if config.args.loci_file is not None:
        fivepseq_out.write_df_to_file(
            fivepseq_counts.get_pauses_from_loci(config.args.loci_file, config.span_size,
                                                 20),
            fivepseq_out.LOCI_PAUSES_FILE)

    # TODO save FivePseqCounts object as pickle
    return fivepseq_counts


def visualize(fivepseq_counts):
    # scatterplots.plot_frames_term(fivepseq_counts, FivePSeqCounts.TERM,
    # os.path.join(config.out_dir, "meta_count_frames.pdf"))
    # scatterplots.plot_frames_term(fivepseq_counts, FivePSeqCounts.START,
    # os.path.join(config.out_dir, "meta_count_frames_start.pdf"))

    pass


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
    fivepseq_counts = generate_and_store_fivepseq_counts()
    visualize(fivepseq_counts)

    # unique_sequences = fivepseq_counts.get_unique_sequences(FivePSeqCounts.TERM, 3, 0)

    # wrap-up
    elapsed_time = time.clock() - start_time
    config.logger.info("SUCCESS! Fivepseq finished in\t%s\tseconds. The report files maybe accessed at:\t\t%s "
                       % (elapsed_time, config.out_dir))

    # TODO handle KeyboardInterrupt's (try, wait, except, pass?)


if __name__ == '__main__':
    main()
