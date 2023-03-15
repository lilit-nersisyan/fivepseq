from __future__ import print_function

"""This module contains classes that are specialized for opening specific types of files.
They also retrieve and store the file properties (e.g. compression, extension) for further reference.
"""
import logging
import os

# PORT: pathlib2 is for python version 2.7, use pathlib in version 3  <>
import dill
import pathlib2
import plastid
import pysam
from preconditions import preconditions

from fivepseq import config
from fivepseq.logic.structures.alignment import Alignment
from fivepseq.logic.structures.annotation import Annotation
from fivepseq.logic.structures.genome import Genome

COMPRESSION_GZ = ".gz"
COMPRESSION_BZ = ".bz"
COMPRESSION_None = ""

logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

class TopReader:
    valid_extensions = None
    file_basename = None
    extension = None
    compression = None
    file = None
    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

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
            raise IOError("Could not instantiate %s. Reason: %s." % (self.__class__.__name__, str(e)))
        except Exception as ex:
            raise Exception("Could not instantiate %s. Reason: %s." % (self.__class__.__name__, str(ex)))

        # check for valid file extension (also store the extension/compression as instance attributes)
        if get_file_extension(file_path) in self.valid_extensions:
            self.compression = get_file_compression(file_path)
            self.extension = get_file_extension(file_path)
            self.file_basename = get_base_file_name(file_path)
        else:
            raise Exception(
                "Could not instantiate %s. Wrong extension of file %s" % (self.__class__.__name__, file_path))

        # set the file as an instance attribute
        self.file = os.path.abspath(file_path)
        self.logger.debug("Initialized %s with file path %s" % (self.__class__.__name__, file_path))

    def __str__(self):
        reader_string = ""
        reader_string += "\tReader: %s\n" % self.__class__.__name__
        reader_string += "\tFile: %s\n" % self.file.name
        reader_string += "\tFile basename: %s\n" % self.file_basename
        reader_string += "\tFile extension: %s\n" % self.extension
        reader_string += "\tFile compression: %s\n" % (
            "None" if self.compression == COMPRESSION_None else self.compression)
        return reader_string


class BamReader(TopReader):
    # TODO should we allow for sam files?
    EXTENSION_BAM = "bam"
    EXTENSION_SAM = "sam"
    three_prime = False
    valid_extensions = [EXTENSION_BAM, EXTENSION_SAM]

    alignment = None

    def __init__(self, file_path, three_prime = False):
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
        self.three_prime = three_prime
        self.create_bam_array()
        # TODO check if index file exists


    def create_bam_array(self):
        """
        Reads the alignment file and converts it to plastid-specific array of
        :return:
        """
        # REQ this will work only if the index file is in the same directory as the bam file, and is named in a standard manner : https://pysam.readthedocs.io/en/latest/api.html
        alignment_file = pysam.AlignmentFile(self.file)
        self.alignment = Alignment(alignment_file, self.file, self.three_prime)
        # TODO the gene set file should be converted to a filter and then used to filter the transcript assembly
        # TODO (or should it be added to the mapping function?)


class FastaReader(TopReader):
    EXTENSION_FA = "fa"
    EXTENSION_FNA = "fna"
    EXTENSION_FASTA = "fasta"
    valid_extensions = [EXTENSION_FA, EXTENSION_FNA, EXTENSION_FASTA]
    genome = None

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
        self.logger.info("Reading in genome file...")
        self.genome = Genome(file_path)


class AnnotationReader(TopReader):
    EXTENSION_GTF = "gtf"
    EXTENSION_GFF = "gff"
    EXTENSION_GFF3 = "gff3"
    valid_extensions = [EXTENSION_GTF, EXTENSION_GFF, EXTENSION_GFF3]
    annotation = None
    CODING_TYPES = ["mRNA", "coding", "protein_coding"]
    CODING = "coding"
    PROTEIN_CODING = "protein_coding"
    MRNA = "mRNA"

    def __init__(self, file_path, transcript_type = "mRNA", break_after=None):
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
            transcript_assembly = self.create_transcript_assembly(biotype = transcript_type, break_after=break_after)

        except Exception as e:
            error_message = "Problem generating transcript assembly from annotation file %s. Reason:%s" % (
                self.file, str(e))
            self.logger.error(error_message)
            raise Exception(error_message)

        self.annotation = Annotation(transcript_assembly, file_path, transcript_type = transcript_type)

    def create_transcript_assembly(self, biotype, break_after=None):
        """
        Uses the plastid transcript assembly generator to retrieve transcripts from the annotation file.
        Stores the transcripts in a list, to be accessed repeatedly later.

        :param biotype:
        :param break_after: int, for testing purposes only: read in a limited amount of transcripts and break
        :return: list of transcripts read from the annotation file
        """

        if biotype is None:
            print("No biotype passed to AnnotationReader. Sticking to default one: %s" % self.CODING)
            biotype = self.CODING
        # attempt to load assembly from pickle path if it exists
        if not config.args.ignore_cache:
            transcript_assembly = self.load_assembly_from_pickle(biotype, break_after)
            if transcript_assembly is not None:
                return transcript_assembly


        # if transcript assembly not loaded proceed to reading it

        self.logger.info("Reading in transcript assembly with transcript type %s" % biotype)

        if self.extension == self.EXTENSION_GTF:
            transcript_assembly_generator = plastid.GTF2_TranscriptAssembler(self.file, return_type=plastid.Transcript)
        else:
            transcript_assembly_generator = plastid.GFF3_TranscriptAssembler(self.file, return_type=plastid.Transcript)

        transcript_assembly = [None] * 1000000
        index = 0
        i = 1
        progress_bar = ""

        biotypes_file = None
        if config.cache_dir is not None:
            biotypes_file = os.path.join(config.cache_dir, os.path.basename(self.file) + "_unique-types.txt")
        unique_biotypes = []

        for transcript in transcript_assembly_generator:
            if break_after is not None:
                if i == break_after:
                    break

            if i % 1000 == 0:
                self.logger.info("\r>>Transcript count: %d\tValid transcripts: %d\t%s" % (i, index, progress_bar), )

            # filter for biotype
            valid_type = False
            if 'biotype' in transcript.attr:
                transcript_biotype = transcript.attr.get('biotype')
            else:
                transcript_biotype = transcript.attr.get('type')


            if self.biotype_valid(biotype, transcript_biotype):
                if transcript.cds_start is not None:
                    valid_type = True
                    if index < len(transcript_assembly):
                        transcript_assembly[index] = transcript
                    else:
                        transcript_assembly.append(transcript)
                    index += 1

                if transcript_biotype not in unique_biotypes:
                    unique_biotypes.append(transcript_biotype)
                    print(transcript_biotype + "\tValid: %s" % valid_type)

            i += 1

        if index == 0:
            error_message = "No transcripts with biotype %s were present" % biotype
            self.logger.error(error_message)
            raise Exception(error_message)

        if biotypes_file is not None:
            with open(biotypes_file, 'w') as file:
                file.write("\n".join(unique_biotypes))

        if index < len(transcript_assembly):
            transcript_assembly = transcript_assembly[:index]


        self.logger.debug(
            "Read %d transcripts, of which %d were of type %s" % (i, index, biotype))

        self.logger.debug("Transcript assembly read to memory, with %d valid transcripts" % index)

        self.save_assembly_to_pickle(transcript_assembly, biotype = biotype, break_after=break_after)

        return transcript_assembly

    def biotype_valid(self, biotype, transcript_biotype):
        if biotype == self.CODING or biotype == self.MRNA or biotype == self.PROTEIN_CODING:
            for type_token in self.CODING_TYPES:
                if type_token in transcript_biotype:
                    return True
        else:
            if transcript_biotype == biotype:
                return True
        return False

    def load_assembly_from_pickle(self, biotype = CODING, break_after = None):
        """
        Loads the transcript assembly from a pickle path if it exists


        :param biotype:
        :param break_after:
        :return: transcript assembly or None if no pickle path exists
        """

        transcript_assembly = None

        if config.cache_dir is not None:
            pickle_path = self.get_assembly_pickle_path(biotype, break_after)

            if pickle_path is not None and os.path.exists(pickle_path):
                self.logger.debug("Loading transcript assembly from existing pickle path %s" % pickle_path)

                try:
                    transcript_assembly = dill.load((open(pickle_path, "rb")))
                    self.logger.debug("Successfully loaded transcript assembly")

                except Exception as e:
                    warning_message = "Problem loading transcript assembly from pickle path %s. Reason: %s" % (
                        pickle_path, str(e))
                    self.logger.warning(warning_message)

        return transcript_assembly

    def get_assembly_pickle_path(self, biotype = CODING, break_after = None):
        """
        Returns the pickle path to an assembly with specified paramaters if such a path exists.

        :param biotype:
        :param break_after:
        :return: the path to existing or non-existing pickle path
        """

        if break_after is None:
            pickle_path = os.path.join(config.cache_dir, os.path.basename(self.file) + "_" + biotype + ".sav")
        else:
            pickle_path = os.path.join(config.cache_dir, os.path.basename(self.file) + "_" + biotype + "_break-" +
                                       str(break_after) + ".sav")

        return pickle_path

    def save_assembly_to_pickle(self, transcript_assembly, biotype = CODING, break_after = None):
        """
        Saves the generated transcript assembly to a pickle path with specified parameters.

        :param transcript_assembly:
        :return:
        """
        pickle_path = self.get_assembly_pickle_path(biotype, break_after)

        if config.cache_dir is not None:
            with open(pickle_path, "wb") as dill_file:
                dill.dump(transcript_assembly, dill_file)
            self.logger.debug("Dumped transcript assembly to file %s" % pickle_path)


@preconditions(lambda file_path: isinstance(file_path, str))
def check_file_path(file_path):
    """
    Checks if the file exists. If it does checks if it is not a directory. Raises IOError otherwise.

    :return:
    """
    if not isinstance(file_path, str):
        raise TypeError("The file path argument is of type %s, expected string" % type(file_path))
    file = pathlib2.Path(file_path)
    if not file.is_file():
        if file.is_dir():
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
