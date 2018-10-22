"""This module contains classes that are specialized for opening specific types of files.
They also retrieve and store the file properties (e.g. compression, extension) for further reference.
"""
import os
# PORT: pathlib2 is for python version 2.7, use pathlib in version 3  <>
import sys

import pathlib2
import plastid
import pysam
from preconditions import preconditions

import fivepseq.config
from fivepseq.logic.structures.alignment import Alignment
from fivepseq.logic.structures.annotation import Annotation

COMPRESSION_GZ = ".gz"
COMPRESSION_BZ = ".bz"
COMPRESSION_None = ""


class TopReader:
    valid_extensions = None
    file_basename = None
    extension = None
    compression = None
    file = None

    @preconditions(lambda file_path: isinstance(file_path, str))
    def __init__(self, file_path):
        """
        Initializes a TopReader with the given file path.
        Checks validity of the file. Raises IOError if the file does not exist or is a directory.
        Checks the file type through its extension. Raises Exception for invalid extensions.
        Stores file compression and extension as attributes.
        :param file_path
        """

        # check if file exists
        try:
            check_file_path(file_path)
        except IOError as e:
            raise IOError("Could not instantiate %s. Reason: %s." % (self.__class__.__name__, e.message))
        except Exception as ex:
            raise Exception("Could not instantiate %s. Reason: %s." % (self.__class__.__name__, ex.message))

        # check for valid file extension (also store the extension/compression as instance attributes)
        if get_file_extension(file_path) in self.valid_extensions:
            self.compression = get_file_compression(file_path)
            self.extension = get_file_extension(file_path)
            self.file_basename = get_base_file_name(file_path)
        else:
            raise Exception("Could not instantiate %s. Wrong extension of file %s" % (self.__class__.__name__, file_path))

        # set the file as an instance attribute
        self.file = os.path.abspath(file_path)
        fivepseq.config.logger.debug("Initialized %s with file path %s" % (self.__class__.__name__, file_path))

    def __str__(self):
        reader_string = ""
        reader_string += "\tReader: %s\n" % self.__class__.__name__
        reader_string += "\tFile: %s\n" % self.file.name
        reader_string += "\tFile basename: %s\n" % self.file_basename
        reader_string += "\tFile extension: %s\n" % self.extension
        reader_string += "\tFile compression: %s\n" % ("None" if self.compression == COMPRESSION_None else self.compression )
        return reader_string


class BamReader(TopReader):
    # TODO should we allow for sam files?
    EXTENSION_BAM = "bam"
    EXTENSION_SAM = "sam"
    valid_extensions = [EXTENSION_BAM, EXTENSION_SAM]

    alignment = None

    def __init__(self, file_path):
        """
        Initializes a BamReader with the given file path.
        Checks the validity of the file. Raises IOError if the file does not exist or is a directory.

        Checks the file type through its extension. Raises Exception for invalid extensions.
        Valid extensions/compressions are:
            .bam    .bam.gz .bam.bz
            .sam

        Stores file compression and extension as attributes.
        :param file_path:
        """
        TopReader.__init__(self, file_path)
        self.create_bam_array()
        # TODO check if index file exists

    def create_bam_array(self):
        """
        Reads the alignment file and converts it to plastid-specific array of
        :return:
        """
        # REQ this will work only if the index file is in the same directory as the bam file, and is named in a standard manner : https://pysam.readthedocs.io/en/latest/api.html
        alignment_file = pysam.AlignmentFile(self.file)
        self.alignment = Alignment(alignment_file)
        # TODO the gene set file should be converted to a filter and then used to filter the transcript assembly
        # TODO (or should it be added to the mapping function?)


class FastaReader(TopReader):
    EXTENSION_FA = "fa"
    EXTENSION_FASTA = "fasta"
    valid_extensions = [EXTENSION_FA, EXTENSION_FASTA]

    def __init__(self, file_path):
        """
        Initializes a FastaReader with the given file path.
        Checks the validity of the file. Raises IOError if the file does not exist or is a directory.

        Checks the file type through its extension. Raises Exception for invalid extensions.
        Valid extensions/compressions are:
            .fa     .fa.gz      .fa.bz
            .fasta  .fasta.gz   .fasta.bz

        Stores file compression and extension as attributes.
        :param file_path:
        """
        TopReader.__init__(self, file_path)

class AnnotationReader(TopReader):
    EXTENSION_GTF = "gtf"
    EXTENSION_GFF = "gff"
    EXTENSION_GFF3 = "gff3"
    valid_extensions = [EXTENSION_GTF, EXTENSION_GFF, EXTENSION_GFF3]
    annotation = None

    def __init__(self, file_path):
        """
        Initializes an AnnotationReader with the given file path.
        Checks the validity of the file. Raises IOError if the file does not exist or is a directory.

        Checks the file type through its extension. Raises Exception for invalid extensions.
        Valid extensions/compressions are:
            .gtf .gtf.gz   .gtf.bz
            .gff .gff.gz   .gff.bz
            .gff3 .gff3.gz   .gff3.bz

        Stores file compression and extension as attributes.

        Creates a plastid transcript assembly and initiates an Annotation object.
        Raises Exception  when transcript assembly generation is not successful.
        :param file_path:
        """

        TopReader.__init__(self, file_path)

        try:
            transcript_assembly= self.create_transcript_assembly()
        except Exception as e:
            error_message = "Problem generating transcript assembly from annotation file %s. Reason:%s" % (self.file, e.message)
            fivepseq.config.logger.error(error_message)
            raise Exception(error_message)

        self.annotation = Annotation(transcript_assembly, file_path)



    def create_transcript_assembly(self):
        """
        Uses the plastid transcript assembly generator to retrieve transcripts from the annotation file.
        Stores the transcripts in a list, to be accessed repeatedly later.
        :return:
        """
        # TODO dill.dump and load
        fivepseq.config.logger.debug("Reading in transcript assembly...")
        if self.extension == self.EXTENSION_GTF:
            transcript_assembly_generator = plastid.GTF2_TranscriptAssembler(self.file, return_type=plastid.Transcript)
        else:
            transcript_assembly_generator = plastid.GFF3_TranscriptAssembler(self.file, return_type=plastid.Transcript)
        transcript_assembly = []
        i = 1
        progress_bar = ""
        for transcript in transcript_assembly_generator:
            if i % 100 == 0:
                progress_bar += "#"
                print "\r>>Transcript count: %d\t%s" % (i, progress_bar),
            transcript_assembly.append(transcript)
            i += 1
        return transcript_assembly

@preconditions(lambda file_path: isinstance(file_path, str))
def check_file_path(file_path):
    """
    Checks if the file exists. If it does checks if it is not a directory. Raises IOError otherwise.

    :return:
    """
    if not isinstance(file_path, str):
        raise TypeError("The file path argument is of type %s, expected string" % type(file_path))
    alignment_file = pathlib2.Path(file_path)
    if not alignment_file.is_file():
        if alignment_file.is_dir():
            raise IOError("The provided alignment file path %s is a directory" % file_path)
        else:
            raise IOError("Could not open %s: no such file or directory" % file_path)


@preconditions(lambda file_path: isinstance(file_path, str))
def get_base_file_name(file_path):
    """
    Returns the base file name without extension and compression.
    :param file_path: file path
    :return: base file name without extension and compression
    """
    compression = get_file_compression(file_path)
    extension = get_file_extension(file_path)
    if os.path.isfile(file_path):
        file_name = os.path.basename(file_path.replace(compression, "").replace("." + extension, ""))
    else:
        raise IOError("The file %s does not exist or is not a file." % file_path)
    return file_name


@preconditions(lambda file_path: isinstance(file_path, str))
def get_file_extension(file_path):
    """
    Returns file extension excluding compression.
        - If the file is uncompressed, returns the string followed by the last dot.
        - If the file is compressed, returns the string followed by the last dot coming before ".gz" or ".bz" compressions.
    :param file_path: the path to the file
    :return: file extension as string
    """
    compression = get_file_compression(file_path)
    file_path_raw = file_path.replace(compression, "")
    if len(file_path_raw.split(".")) < 2:
        return ""
    else:
        return file_path_raw.split(".")[-1]


@preconditions(lambda file_path: isinstance(file_path, str))
def get_file_compression(file_path):
    """
    Returns the compression of the file if is one of ".bz" or ".bz".
    Returns an empty string otherwise
    :param file_path: file path
    :return: compression string, either ".bz", ".gz" or ""
    """
    if file_path.endswith(".gz"):
        return COMPRESSION_GZ
    elif file_path.endswith(".bz"):
        return COMPRESSION_BZ
    else:
        return COMPRESSION_None
