from __future__ import print_function

import fivepseq

"""
The fivepseq entry point.
"""
# PORT
# argparse was written for argparse in Python 3.
# A few details are different in 2.x, especially some exception messages, which were improved in 3.x.
import argparse
import glob
import logging
import os
import sys
import time

from fivepseq import config
from fivepseq.logic.algorithms.general_pipelines.count_pipeline import CountPipeline
from fivepseq.logic.algorithms.general_pipelines.viz_pipeline import VizPipeline
from fivepseq.logic.structures.annotation import Annotation
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts, CountManager
from fivepseq.util.formatting import pad_spaces
from fivepseq.util.readers import BamReader, AnnotationReader, FastaReader
from fivepseq.util.writers import FivePSeqOut

FIVEPSEQ_COUNTS_DIR = "fivepseq_counts"
FIVEPSEQ_PLOTS_DIR = "fivepseq_plots"
FIVEPSEQ_LOG_DIR = "log"


class FivepseqArguments:
    """
    This class is for reading, parsing, storing and printing command line arguments.
    """
    logger = None

    def __init__(self):
        """
        Initializes FivepseqArguments with command line arguments.
        Stores those in the <args> attribute.
        """

        ###############################
        #       fivepseq common args
        ###############################

        parser = argparse.ArgumentParser(prog="fivepseq",
                                         description="Reports the distribution of 5\' endpoints in RNA-seq reads from alignment files.")
        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        additional = parser.add_argument_group('additional arguments')
        advanced = parser.add_argument_group('advanced arguments')

        optional.add_argument('--version', action='version', version='%(prog)s' + fivepseq.__version__)
        # TODO accept compressed files

        ##########################
        #       required arguments
        ##########################

        required.add_argument("-b", "-bam",
                              help="the full path one or many bam/sam files (many files should be provided with a pattern, within quotation marks, e.g. \"bam/*.bam\")",
                              type=str,
                              required=True,
                              metavar="[full path to BAMs with pattern in quotation marks, e.g. \"bam/*.bam\"]")

        required.add_argument("-g", "-genome",
                              help="the full path to the fa/fasta file",
                              type=str,
                              required=True,
                              metavar="FASTA")

        required.add_argument("-a", "-annotation",
                              help="the full path to the gtf/gff/gff3 file",
                              type=str,
                              required=True,
                              metavar="GFF")

        optional.add_argument("-o", "-outdir",
                              help="the output directory",
                              type=str,
                              required=False,
                              default="fivepseq_out",
                              metavar="output_directory")

        optional.add_argument("-t", "-title",
                              help="title of the reports",
                              type=str,
                              required=False,
                              metavar="title")

        optional.add_argument("--span",
                              help="the number of bases to span around a genomic position",
                              type=int,
                              required=False,
                              choices=range(0, 1000),
                              default=100)

        optional.add_argument('--conflicts',
                              help="file conflict mode:\n"
                                   "add (default) - only adds missing files to existing directory\n"
                                   "overwrite - overwrites all the files in the existing (in case) output directory\n"
                                   "alt_dir - uses alternative directory by appending '+' suffix to existing (in case) output directory",
                              type=str,
                              choices=['add', 'overwrite', 'alt_dir'],
                              required=False,
                              default="add")

        additional.add_argument("-gf", "-genefilter",
                                help="file containing gene names of interest",
                                type=str,
                                required=False)

        additional.add_argument("-gs", "-genesets",
                                help="a file with gene-geneset mapping",
                                type=str,
                                required=False)

        optional.add_argument("--ignore-cache",
                              help="read in transcript assembly and rewrite existing pickle paths",
                              action="store_true",
                              default=False,
                              required=False)

        optional.add_argument("-transcript-type",
                              help="choose transcript type (attribute biotype) to analyze",
                              type=str,
                              default=Annotation.default_filter,
                              required=False)

        optional.add_argument("--no-mask",
                              help="Turn off masking transcript boundaries for codon-relative counts",
                              action="store_true",
                              default=False,
                              required=False)

        optional.add_argument("-subset",
                              help="Use a given size subset of transcripts",
                              required=False,
                              type=int
                              )

        advanced.add_argument('--thrp', '--three-prime',
                              help="Should the 3' ends used as datapoints instead of 5' default?",
                              action="store_true",
                              default=False,
                              required=False)

        advanced.add_argument('-ms', "-codon-mask-size",
                              help="a number {0, >3} specifying how many positions to mask from transcript start and end when counting codon-relative counts",
                              type=int,
                              required=False,
                              default=20,
                              choices=range(3, 50),
                              metavar="[3:50)")

        advanced.add_argument("-dipeptide-pos",
                              help="counts in which position {-27:6} from A site should be ordered to output top stalled dipeptides",
                              type=int,
                              required=False,
                              default=-14,
                              choices=range(-27, 6),
                              metavar="[-27:6)"
                              )

        advanced.add_argument("-tripeptide-pos",
                              help="counts in which position {-24:9} from A site should be ordered to output top stalled tripeptides",
                              type=int,
                              required=False,
                              default=-11,
                              choices=range(-24, 9),
                              metavar="[-24:9)"
                              )

        advanced.add_argument("-triangle_threshold",
                              help="count threshold for points in the triangle plot (points with lower counts will not be plotted)",
                              type=int,
                              required=False,
                              default=10,
                              choices=range(0, 10000),
                              metavar="0:10000"
                              )

        advanced.add_argument("--loci-file",
                              help="coordinates of loci relative to which mapping positions are to be plotted",
                              required=False,
                              type=str)

        advanced.add_argument('--ds', '--downsample',
                              help="constant threshold for down-sampling: points exceeding this value will be downsampled"
                                   "to be equal to it",
                              type=float,
                              required=False)

        advanced.add_argument('-op', '-outlier_probability',
                              help="probablity threshold for poisson distribution: points less than this value will be downsampled"
                                   "as outliers",
                              type=float,
                              required=False,
                              default=0)

        advanced.add_argument("--log",
                              help="set logging level",
                              choices=["INFO", "DEBUG", "info", "debug"],
                              default="DEBUG",
                              required=False)

        advanced.add_argument("-queue",
                              # help="perform rq-score analysis with the specified size",
                              type=int,
                              # default=30,
                              required=False)

        advanced.add_argument("-cspan",
                              # help="perform rq-score analysis with the specified size",
                              type=int,
                              default=-30,
                              required=False)

        advanced.add_argument("--oof",
                              help="compute out-of-frame counts",
                              action="store_true",
                              default=False,
                              required=False)

        advanced.add_argument("-fpat", "-filter_pattern",
                              help="reads proceeding this sequence pattern will be skipped",
                              type=str,
                              required=False)

        advanced.add_argument("--aa-islands",
                              help="amino acid line charts are drawn for island regions for each amino acid, in a cspan range (default -30:+6)",
                              action="store_true",
                              default=False,
                              required=False)

        #        parser.add_argument("-tf", "-transcript_filter",
        #                            help="Name of filter to apply on transcripts",
        #                            type=str,
        #                            choices=[VizPipeline.FILTER_TOP_POPULATED,
        #                                     VizPipeline.FILTER_CANONICAL_TRANSCRIPTS,
        #                                     Annotation.FORWARD_STRAND_FILTER,
        #                                     Annotation.REVERSE_STRAND_FILTER],
        #                            required=False)

        ##########################
        #       config setup
        ##########################

        config.args = parser.parse_args()
        #        if config.args.command == 'count':
        #            config.bam = os.path.abspath(config.args.b)
        if not hasattr(config.args, "command"):
            config.args.command = "count_and_plot"

        if not hasattr(config.args, "transcript_type") or config.args.transcript_type is None:
            config.args.transcript_type = "protein_coding"

        if hasattr(config.args, "no_mask") and config.args.no_mask:
            config.args.codon_mask_size = 0
        elif not hasattr(config.args, "codon_mask_size"):
            config.args.codon_mask_size = 20

        # FIXME this imposes necessity to know the span size with which fivepseq was run for further visualization
        config.span_size = config.args.span
        config.genome = os.path.abspath(config.args.g)
        config.annot = os.path.abspath(config.args.a)
        config.out_dir = os.path.abspath(config.args.o)

    def print_args(self):
        """
        Prints the command line arguments in a formatted manner.
        :return:
        """
        # NOTE The files are not checked for existence at this stage.
        # NOTE The absolute_path function simply merges the relative paths with the working directory.

        if config.args.command == 'count':
            print("Call to fivepseq count with arguments:\n")
            print("%s%s" % (pad_spaces("\tAlignment file:"), os.path.abspath(config.args.b)))
            print("%s%s" % (pad_spaces("\tGenome file:"), os.path.abspath(config.args.g)))
            print("%s%s" % (pad_spaces("\tAnnotation file:"), os.path.abspath(config.args.a)))

            print("%s%s" % (pad_spaces("\tSpan size:"), config.span_size))
            if config.args.gf is not None:
                print("%s%s" % (pad_spaces("\tGenefilter file:"), os.path.abspath(config.args.gf)))
            if config.args.gs is not None:
                print("%s%s" % (pad_spaces("\tFunctional gene sets file:"), os.path.abspath(config.args.gs)))
            print("%s%s" % (pad_spaces("\tOutput directory:"), os.path.abspath(config.args.o)))
            print("%s%s" % (pad_spaces("\tConflict handling:"), config.args.conflicts))

            # advanced arguments
            print("%s%s" % (pad_spaces("\tTranscript boundary masking for codon counts:"), config.args.no_mask))
            if not config.args.no_mask:
                print("%s%s" % (
                    pad_spaces("\tTranscript boundary mask size for codon counts:"), config.args.codon_mask_size))

            print("%s" % (pad_spaces(
                "\tTripeptides will be sorted at position %d from the A site:" % config.args.tripeptide_pos)))

            print("%s" % (pad_spaces(
                "\tDipeptides will be sorted at position %d from the A site:" % config.args.dipeptide_pos)))

            if hasattr(config.args, "triangle_threshold") and config.args.triangle_threshold is not None:
                print("%s" % (pad_spaces(
                    "\tGenes with more than %s counts will be included in the triangle plot:" %
                    config.args.triangle_threshold)))


        elif config.args.command == 'plot':
            config.args.count_folders = []
            print("Call to fivepseq plot with arguments:\n")
            if config.args.sd is not None:
                print("%s%s" % (pad_spaces("\tInput dir:"), os.path.abspath(config.args.sd)))
                config.args.count_folders.append(os.path.abspath(config.args.sd))
            else:
                if len(glob.glob(config.args.md)) == 0:
                    err_msg = "Provided input %s does not contain proper fivepseq directories" % config.args.md
                    raise Exception(err_msg)
                print("%s" % (pad_spaces("\tInput directories:")))

                # TODO check for proper patterned input
                for f in glob.glob(config.args.md):
                    config.args.count_folders.append(f)
                    print("%s" % pad_spaces("\t%s" % f))
            print("%s%s" % (pad_spaces("\tGenome file:"), os.path.abspath(config.args.g)))
            print("%s%s" % (pad_spaces("\tAnnotation file:"), os.path.abspath(config.args.a)))
            if config.args.o is not None:
                print("%s%s" % (pad_spaces("\tOutput directory:"), os.path.abspath(config.args.o)))

            if config.args.t is not None:
                print("%s%s" % (pad_spaces("\tOutput file title:"), config.args.t + ".html"))

                # NOTE The files are not checked for existence at this stage.
                # NOTE The absolute_path function simply merges the relative paths with the working directory.

        elif config.args.command == 'count_and_plot':
            print("Call to fivepseq count_and_plot with arguments:\n")
            print("%s%s" % (pad_spaces("\tAlignment file:"), os.path.abspath(config.args.b)))
            print("%s%s" % (pad_spaces("\tGenome file:"), os.path.abspath(config.args.g)))
            print("%s%s" % (pad_spaces("\tAnnotation file:"), os.path.abspath(config.args.a)))

            print("%s%s" % (pad_spaces("\tSpan size:"), config.span_size))
            if config.args.gf is not None:
                print("%s%s" % (pad_spaces("\tGenefilter file:"), os.path.abspath(config.args.gf)))
            if config.args.gs is not None:
                print("%s%s" % (pad_spaces("\tGenesets file:"), os.path.abspath(config.args.gs)))
            print("%s%s" % (pad_spaces("\tOutput directory:"), os.path.abspath(config.args.o)))
            print("%s%s" % (pad_spaces("\tConflict handling:"), config.args.conflicts))

            print("%s%s" % (pad_spaces("\tLogging level:"), config.args.log.upper()))

            if config.args.t is not None:
                print("%s%s" % (pad_spaces("\tOutput file title:"), config.args.t + ".html"))

            if hasattr(config.args, "subset") and config.args.subset is not None:
                print("Number of transcripts to subset: %s" % config.args.subset)
            # optional arguments
            if hasattr(config.args, "transcript_type"):
                print("%s%s" % (pad_spaces("\tTranscript type:"), config.args.transcript_type))
            # advanced arguments
            if hasattr(config.args, "thrp") and config.args.thrp:
                print("%s" % (pad_spaces("\t3' prime mode activated!")))

            if config.args.no_mask:
                print("%s%s" % (pad_spaces("\tTranscript boundary masking for codon counts:"), "OFF"))
            else:
                print("%s%s" % (
                    pad_spaces("\tTranscript boundary mask size for codon counts:"), config.args.codon_mask_size))

            if hasattr(config.args, "tripeptide_pos"):
                print("%s%s" % (pad_spaces("\tSort tripeptides from A site at:"), config.args.tripeptide_pos))

            if hasattr(config.args, "dipeptide_pos"):
                print("%s%s" % (pad_spaces("\tSort dipeptides from A site at:"), config.args.dipeptide_pos))

            if hasattr(config.args, "triangle_threshold") and config.args.triangle_threshold is not None:
                print("%s%s" % (pad_spaces("\tCounts threshold to include genes in the triangle plot:"),
                                config.args.triangle_threshold))

            if hasattr(config.args, "queue") and config.args.queue is not None:
                print("%s%s" % (pad_spaces("\tRibosome queue score will be computed with size:"),
                                config.args.queue))

            if hasattr(config.args, "cspan"):
                print("%s%s" % (pad_spaces("\tspan size for codon counts:"), config.args.cspan))

            if config.args.loci_file is not None:
                print("%s%s" % (pad_spaces("\tLoci file:"), config.args.loci_file))

            if config.args.fpat is not None:
                print("%s%s" % (pad_spaces("\tFilter pattern:"), config.args.fpat))

            if hasattr(config.args, "аа_islands") and config.args.аа_islands:
                print("%s" % (pad_spaces("\tAmino acid linecharts for islands will be plotted!")))

        #   if config.args.fivepseq_pickle is not None:
        #       print "%s%s" % (pad_spaces("\tFivepseq pickle path specified:"), config.args.fivepseq_pickle)


def setup_output_dir():
    """
    This method attempts to set the output directory for fivepseq count reports.

    If the user-specified output_dir does not exist and the parent directory does, then the output will be stored there.

    If the user-specified output_dir does not exist and the parent directory also, then an IOError is raised.

    If the user-specified output_dir already exists, then the file conflicts are resolved based on the option specified:
        - add (default): only the count files not in the directory are generated and added
        - overwrite:  attempts to remove and replace all the content of the directory. Exception is raised in case of failure.
        - alt_dir: appends "+" signs to the  end of the directory name recursively, until no such directory exists. Creates a new directory and stores the output there.

    :return:
    """

    output_dir = config.out_dir
    if os.path.exists(output_dir):
        if config.args.conflicts == config.ADD_FILES:
            print("\nWARNING:\tThe output directory %s already exists.\n"
                  "\tOnly missing files will be added." % output_dir)

        if config.args.conflicts == config.OVERWRITE:
            print("\nWARNING:\tThe output directory %s already exists.\n"
                  "\tOnly missing files will be added." % output_dir)

        elif config.args.conflicts == config.ALT_DIR:
            while os.path.exists(output_dir):
                output_dir += "+"
            print("\nWARNING:\tOutput directory %s already exists.\n\tFiles will be stored in %s" % (
                config.out_dir, output_dir))
            config.out_dir = output_dir
            make_dir(output_dir)

    else:
        make_dir(output_dir)

    # FIXME probably the cache directory should reside somewhere in the fivepseq package folder
    # TODO add checks to make sure the parent directory exists and has write permissions
    # TODO (maybe should not be a problem when it is within the fivepseq package, or not?)
    config.cache_dir = os.path.join(os.path.dirname(config.out_dir), "fivepseq_cache")
    if not os.path.exists(config.cache_dir):
        os.makedirs(config.cache_dir)

    # create log directory
    if not os.path.exists(os.path.join(config.out_dir, FIVEPSEQ_LOG_DIR)):
        try:
            os.mkdir(os.path.join(config.out_dir, FIVEPSEQ_LOG_DIR))
        except Exception as e:
            err_msg = "Problem creating log directory: %s" % str(e)
            raise Exception(err_msg)

    # create fivepseq_counts and fivepseq_plots directories for count_and_plot case
    if config.args.command == 'count_and_plot':
        if not os.path.exists(os.path.join(config.out_dir, FIVEPSEQ_COUNTS_DIR)):
            try:
                os.mkdir(os.path.join(config.out_dir, FIVEPSEQ_COUNTS_DIR))
            except Exception as e:
                err_msg = "Problem making directory for counts: %s" % str(e)
                logging.getLogger(config.FIVEPSEQ_LOGGER).error(err_msg)
                raise Exception(err_msg)

        if not os.path.exists(os.path.join(config.out_dir, FIVEPSEQ_PLOTS_DIR)):
            try:
                os.mkdir(os.path.join(config.out_dir, FIVEPSEQ_PLOTS_DIR))
            except Exception as e:
                err_msg = "Problem making directory for counts: %s" % str(e)
                logging.getLogger(config.FIVEPSEQ_LOGGER).error(err_msg)
                raise Exception(err_msg)


def make_dir(dir):
    try:
        os.makedirs(dir)
        print("\n\tSETUP:\tSuccess")
        print("\tSETUP:\t%s%s" % (pad_spaces("Output directory:"), dir))
    except Exception as e:
        error_message = "\tSETUP:\tCould not create output directory %s. Reason: %s" % (dir, str(e))
        raise Exception(error_message)


def setup_logger():
    """
    Sets up the logger and the handlers for different logging levels.
    Set the logger level according to the argument --log passed to the argparse. Raises a ValueError in case of an invalid argument.
    Raises an exception if the files for handlers cannot be generated.
    :return:
    """
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                        format='%(levelname)s:%(asctime)s\t [%(filename)s:%(lineno)s - %(funcName)s]\t%(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S')

    fivepseq_logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

    log_file = os.path.join(config.out_dir, FIVEPSEQ_LOG_DIR, "fivepseq.log")
    log_level = getattr(logging, config.args.log.upper(), None)
    if not isinstance(log_level, int):
        raise ValueError('Invalid log level: %s' % config.args.log)

    # INFO level file handler
    log_file_handler = logging.FileHandler(log_file)
    fivepseq_logger.addHandler(log_file_handler)

    if log_level == logging.DEBUG:
        log_file_handler.setLevel(logging.DEBUG)
    else:
        log_file_handler.setLevel(logging.INFO)

    if not os.path.exists(log_file):
        raise Exception("Could not instantiate the logger. Exiting.")

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\tSETUP:\t%s%s"
                                                   % (pad_spaces("Log file:"), log_file))


def generate_and_store_fivepseq_counts(plot=False):
    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)
    logger.info("FIVEPSEQ STARTED")
    logger.info("________________")

    # process bam input

    logger.info("%s" % (pad_spaces("\tInput bam files:")))
    bam_files = []
    for bam in sorted(glob.glob(config.args.b)):
        if bam.endswith(".bam"):
            bam_files.append(bam)
            logger.info("%s" % pad_spaces("\t%s" % bam))
        else:
            logger.info("non bam file found: %s" % bam)

    if len(bam_files) == 0:
        err_msg = "No bam files found at %s" % config.args.b
        logging.getLogger(config.FIVEPSEQ_LOGGER).error(err_msg)
        return None

    # set up annotation

    if hasattr(config.args, "subset") and config.args.subset is not None:
        annotation_reader = AnnotationReader(config.annot, transcript_type=config.args.transcript_type,
                                             break_after=config.args.subset)
    else:
        annotation_reader = AnnotationReader(config.annot, transcript_type=config.args.transcript_type)

    # set the break for human
    annotation = annotation_reader.annotation

    annotation.set_default_span_size(config.span_size)

    if hasattr(config.args, 'gf') and config.args.gf is not None:
        annotation.set_permanent_gene_filter(config.args.gf)

    if hasattr(config.args, 'gs') and config.args.gs is not None:
        annotation.store_gene_sets(config.args.gs)
        if len(annotation.gs_transcriptInd_dict) > 12:
            err_msg = "Too many genesets (%d) provided: should not exceed 12" % len(annotation.gs_transcript_dict)
            logger.error(err_msg)
            raise Exception(err_msg)

    # set up genome

    fasta_reader = FastaReader(config.genome)

    success_values = {}
    # fivepseq_counts_dict = {}
    count_folders = []
    count_folders_gs = []

    for bam in bam_files:
        # set up bam input and output
        if hasattr(config.args, "thrp") and config.args.thrp:
            bam_reader = BamReader(bam, config.args.thrp)
        else:
            bam_reader = BamReader(bam)
        bam_name = os.path.basename(bam)
        if bam_name.endswith(".bam"):
            bam_name = bam_name[0:len(bam_name) - 4]

        bam_out_dir = os.path.join(config.out_dir, FIVEPSEQ_COUNTS_DIR, bam_name)

        if not os.path.exists(bam_out_dir):
            os.mkdir(bam_out_dir)

        main_type_counts = bam_filter_counts(bam_name, bam_reader.alignment,
                                             annotation, fasta_reader.genome,
                                             bam_out_dir, count_folders, success_values,
                                             filter_name=config.args.transcript_type,
                                             loci_file=config.args.loci_file)

        # fivepseq_counts_dict.update({bam: fivepseq_counts})

        # run on genesets
        if annotation.gs_transcriptInd_dict is not None:
            outlier_lower = main_type_counts.get_outlier_lower()
            for gs in annotation.gs_transcriptInd_dict.keys():
                bam_filter_counts(bam_name, bam_reader.alignment, annotation, fasta_reader.genome,
                                  bam_out_dir, count_folders_gs, success_values,
                                  downsample_constant=outlier_lower,
                                  filter_name=gs, filter=annotation.gs_transcriptInd_dict[gs],
                                  is_geneset=True)

    # check if all the files in all directories are in place and store a summary of all runs
    fivepseq_out = FivePSeqOut(config.out_dir,
                               config.OVERWRITE)  # overwrite is for always removing existing summary file

    fivepseq_out.write_dict(success_values, FivePSeqOut.BAM_SUCCESS_SUMMARY)

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\n##################")
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\n# Finished counting successfully! Proceeding to plotting")
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\n##################")

    if plot:
        if len(count_folders) == 0:
            err_msg = "None of the count directories succeeded. Plots will not be generated."
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(err_msg)
        else:
            # set up job title if none is provided
            if not hasattr(config.args, 't') or config.args.t is None:
                config.args.t = os.path.basename(os.path.dirname(config.args.b)) + "_" + os.path.basename(config.args.b)
                config.args.t = config.args.t.replace(".bam", "")
                config.args.t = config.args.t.replace("_*", "")
                config.args.t = config.args.t.replace("*", "")

            # FIXME currently all count folders in the output directory are used for plotting.
            # FIXME this introduces conflicts with pre-existing count files in the folder
            # FIXME the adding of count_folders list shouuld fix for this: need testing
            # config.args.md = str(os.path.join(config.out_dir, FIVEPSEQ_COUNTS_DIR)) + "/*"
            config.args.o = config.out_dir = os.path.join(config.out_dir, FIVEPSEQ_PLOTS_DIR)
            generate_plots(count_folders, annotation.gs_transcriptInd_dict)


def bam_filter_counts(bam_name, alignment, annotation, genome, bam_out_dir,
                      count_folders, success_values, downsample_constant=None,
                      filter_name="protein_coding", is_geneset=False, filter=None, loci_file=None):
    logging.getLogger(config.FIVEPSEQ_LOGGER). \
        info("\n##################\nCounting for sample %s and gene set %s\n##################\n"
             % (bam_name, filter_name))

    filter_out_dir = os.path.join(bam_out_dir, filter_name)
    if not os.path.exists(filter_out_dir):
        os.mkdir(filter_out_dir)

    if filter is not None:
        annotation.apply_geneset_filter(filter_name)

    # combine objects into FivePSeqCounts object
    fivepseq_counts = FivePSeqCounts(alignment, annotation, genome,
                                     # outlier_probability=config.args.op,
                                     config=config,
                                     downsample_constant=downsample_constant,
                                     is_geneset=is_geneset)
    fivepseq_counts.loci_file = loci_file
    fivepseq_counts.fpat = config.args.fpat

    # set up fivepseq out object for this bam
    fivepseq_out = FivePSeqOut(filter_out_dir, config.args.conflicts)

    # run
    fivepseq_pipeline = CountPipeline(fivepseq_counts, fivepseq_out)
    fivepseq_pipeline.run()

    annotation.remove_geneset_filter()

    success = fivepseq_out.sanity_check_for_counts()
    if success:
        success_values.update({bam_name + "_GS_" + filter_name: "SUCCESS"})
        count_folders.append(filter_out_dir)
    else:
        success_values.update({filter_out_dir: "FAILURE"})

    return fivepseq_counts


def generate_plots(count_folders, gs_transcriptInd_dict=None):
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\n#########################")
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\n#  Fivepseq plot called #")
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("\n#########################")

    viz_pipeline = VizPipeline(config.args, count_folders=count_folders, gs_transcriptInd_dict=gs_transcriptInd_dict)
    viz_pipeline.run()

    if hasattr(config.args, 'sd') and config.args.sd is not None:
        print("Single sample %s provided" % config.args.sd)
    else:
        print("Multiple samples in provided")

    pass


def main():
    # argument handling
    fivepseq_arguments = FivepseqArguments()
    fivepseq_arguments.print_args()

    # startup
    setup_output_dir()
    setup_logger()

    start_time = time.time()

    # body

    if (config.args.command == 'count') | (config.args.command == 'count_and_plot'):

        if config.args.command == 'count_and_plot':
            generate_and_store_fivepseq_counts(plot=True)
        else:
            generate_and_store_fivepseq_counts(plot=False)
        elapsed_time = time.time() - start_time

        logging.getLogger(config.FIVEPSEQ_LOGGER).info("SUCCESS! Fivepseq count finished in\t%s\tseconds. "
                                                       "The report files maybe accessed at:\t\t%s "
                                                       % (elapsed_time, config.out_dir))

    elif config.args.command == 'plot':
        plot_logger = logging.getLogger(config.FIVEPSEQ_LOGGER)
        plot_logger.info("Fivepseq plot started")
        generate_plots(config.args.count_folders)
        elapsed_time = time.clock() - start_time
        if config.args.t is not None:
            title = config.args.t
        else:
            if config.args.sd is not None:
                title = os.path.basename(config.args.sd)
            else:
                title = os.path.basename(os.path.dirname(config.args.md)) + "_" + os.path.basename(config.args.md)

        plot_logger.info("SUCCESS! Fivepseq plot finished in\t%s\tseconds. "
                         "The report files maybe accessed at:\t\t%s/%s "
                         % (elapsed_time, config.out_dir, title))

    # TODO handle KeyboardInterrupt's (try, wait, except, pass?)


if __name__ == '__main__':
    main()
