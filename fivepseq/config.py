"""
This module stores singleton objects keeping command line arguments, fivepseq options and loggers.
"""
import logging
import sys

args = None
bam = None
genome = None
annot = None
out_dir = None
cache_dir = None
span_size = None

# default (pre-startup) logger
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(levelname)s:%(asctime)s\t [%(filename)s:%(lineno)s - %(funcName)s]\t%(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S')
logger = logging.getLogger()
