"""
This module keeps properties of annotation files and functions associated with those.
"""
import logging

import plastid
from preconditions import preconditions

from fivepseq import config


class Annotation:
    file_path = None
    geneset_filter = None

    transcript_assembly_dict = {}

    FORWARD_STRAND_FILTER = "fwd"
    REVERSE_STRAND_FILTER = "rev"
    NO_FILTER = "none"

    transcript_filter = NO_FILTER
    span_size = 0
    
    logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)

    @preconditions(lambda file_path: isinstance(file_path, str),
                   lambda transcript_assembly: isinstance(transcript_assembly, list))
    def __init__(self, transcript_assembly, file_path):

        """
        Initiates an Annotation instance with transcript_assembly list.
        Raises a ValueError if transcript assembly is None

        :param transcript_assembly:
        """

        self.file_path = file_path
        if transcript_assembly is None:
            error_message = "transcript_assembly is None. Cannot instantiate an Annotation object."
            self.logger.error(error_message)
            raise ValueError(error_message)
        self.transcript_assembly = transcript_assembly  # the default transcript assembly with no span size and no filter applied
        self.transcript_count = len(transcript_assembly)

        self.transcript_assembly_dict.update({self.NO_FILTER: {
            0: transcript_assembly}})  # the trnascript assembly returned by the reader is set for 0 span size and no filter

        self.logger.debug("Annotation object successfully created from file %s" % file_path)

    def set_default_transcript_filter(self, transcript_filter):
        """
        Sets the default transcript filter to be used with transcript assembly retrieval with no filter specifications.
        """

        self.transcript_filter = transcript_filter

    def set_default_geneset_filter(self, geneset_filter):
        """
        Sets the default gen set filter to be filter the transcript assembly.
        """

        self.geneset_filter= geneset_filter


    def set_default_span_size(self, span_size):
        """
        Sets the default span size to be used with transcript assembly retrieval with no span size specifications.
        """

        self.span_size = span_size

    def get_default_transcript_assembly(self):
        """
        Returns the transcript assembly with default span size and transcript filter.

        :return: list[Transcript] list of transcripts spanned and filtered according to default settings
        """

        return self.get_transcript_assembly(self.span_size, self.transcript_filter)

    def get_transcript_assembly_default_filter(self, span_size = 0):
        """
        Returns the transcript assembly with default filters and specified span size.

        :return: list[Transcript] list of transcripts spanned and filtered according to default settings
        """

        return self.get_transcript_assembly(span_size=0, transcript_filter=self.transcript_filter)

    def get_clean_transcript_assembly(self):
        """
        Returns the transcript assembly as found in the gff file - with no spanning and no transcript filter.

        :return: list[Transcript] list of transcripts spanned and filtered according to default settings
        """

        return self.get_transcript_assembly(span_size=0, transcript_filter=self.NO_FILTER)

    @preconditions(lambda span_size: isinstance(span_size, int))
    def get_transcript_assembly(self, span_size, transcript_filter=None):
        """
        Returns transcripts if a transcript assembly already exists with provided span_size and filter.
        Creates and returns a new transcript assembly otherwise.

        :param span_size:  int span_size - defines the number of nucleotides to span
        :param transcript_filter: filter to be applied on transcripts
        :return: a vector of transcript or an iterable of type Transcript
        """

        if transcript_filter is None:
            transcript_filter = self.NO_FILTER

        if transcript_filter in self.transcript_assembly_dict and \
                span_size in self.transcript_assembly_dict[transcript_filter]:
            
            return self.transcript_assembly_dict[transcript_filter][span_size]

        return self.generate_transcript_assembly(span_size, transcript_filter)

    @preconditions(lambda span_size: isinstance(span_size, int))
    def generate_transcript_assembly(self, span_size, transcript_filter=None):
        """
        Adds spanning regions to the start and the end of the transcript and yields the transcripts one-by-one.
        The original transcript assembly is left intact.

        :param span_size:  int span_size - defines the number of nucleotides to span
        :param transcript_filter: filter to apply on transcripts

        :return: transcript vector
        """

        # error check
        if span_size < 0:
            error_message = "Negative span_size of %d provided: cannot span the transcripts with negative or " \
                            "zero values" % span_size
            self.logger.error(error_message)
            raise ValueError(error_message)
        self.logger.debug("Transcript span size: %d" % span_size)

        # info 
        self.logger.info(
            "Generating transcript assembly under filter: %s and span size: %d"
            % (transcript_filter, span_size))
        
        # generation
        this_transcript_assembly = []
        for transcript in self.transcript_assembly:
            # apply filter
            if transcript_filter == self.FORWARD_STRAND_FILTER:
                if transcript.strand != '+':
                    continue

            elif transcript_filter == self.REVERSE_STRAND_FILTER:
                if transcript.strand != '-':
                    continue

            # span
            spanned_transcript = self.span_transcript(transcript, span_size)

            # append 
            this_transcript_assembly.append(spanned_transcript)
        
        # add assembly to dictionary
        if transcript_filter in self.transcript_assembly_dict:
            self.transcript_assembly_dict[transcript_filter].update({span_size: this_transcript_assembly})
        else:
            self.transcript_assembly_dict.update({transcript_filter: {span_size: this_transcript_assembly}})

        # info
        self.logger.info("Transcript assembly generated")

        return this_transcript_assembly


    @preconditions(lambda span_size: isinstance(span_size, int),
                   lambda span_size: span_size >= 0)
    def span_transcript(self, transcript, span_size):
        """
        Spans the transcript symmetrically by the given span size.
        Spanning may result in going beyond boundaries of current genome, with negative and out of range coordinates,
        which should be handled explicitly in other functions.

        :param transcript:
        :param span_size:
        :return:
        """
        span = transcript.spanning_segment

        new_start = span.start - span_size

        start_span = plastid.GenomicSegment(
            span.chrom, new_start, span.start, span.strand)

        end_span = plastid.GenomicSegment(
            span.chrom, span.end, span.end + span_size, span.strand)

        spanned_transcript = transcript.__deepcopy__(None)
        spanned_transcript.add_segments(start_span, end_span)

        return spanned_transcript