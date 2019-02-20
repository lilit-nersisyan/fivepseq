"""
This module stores singleton objects keeping command line arguments, fivepseq options and loggers.
"""
import logging
import os
import sys

args = None
bam = None
genome = None
annot = None
out_dir = None
cache_dir = None
span_size = None
logger = None

# output conflict modes
ADD_FILES = "add"
OVERWRITE = "overwrite"
ALT_DIR = "alt_dir"

FIVEPSEQ_COUNTS_PICKLE_NAME = "fivepseq.counts.sav"

# default (pre-startup) logger
#logging.basicConfig(stream=sys.stdout, level=logging.INFO,
#                    format='%(levelname)s:%(asctime)s\t [%(filename)s:%(lineno)s - %(funcName)s]\t%(message)s',
#                    datefmt='%m/%d/%Y %I:%M:%S')
#logger = logging.getLogger()

def getFivepseqCountsPicklePath():
    if args.fivepseq_pickle is not None:
        if os.path.exists(args.fivepseq_pickle):
            return args.fivepseq_pickle
        else:
            raise IOError("Specified fivepseq pickle path %s does not exist" % args.fivepseq_pickle)

    pickle_path = os.path.join(out_dir, FIVEPSEQ_COUNTS_PICKLE_NAME)
    if os.path.exists(pickle_path):
        return pickle_path
    if args.conflicts == ALT_DIR:
        while out_dir[-1] == '+':
            prev_dir = out_dir[0:(len(out_dir)-1)]
            if os.path.exists(prev_dir):
                pickle_path = os.path.join(prev_dir, FIVEPSEQ_COUNTS_PICKLE_NAME)
                if os.path.exists(pickle_path):
                    return pickle_path
    return None


def generateFivepseqPicklePath():
    pickle_path = os.path.join(out_dir, FIVEPSEQ_COUNTS_PICKLE_NAME)
    return pickle_path

