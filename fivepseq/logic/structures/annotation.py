"""
This module keeps properties of annotation files and functions associated with those.
"""
import logging

import plastid
from preconditions import preconditions

import config
import fivepseq.config


class Annotation:
    file_path = None
    transcript_assembly = None
    geneset_filter = None
    transcript_count = None
    logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)

    transcript_assembly_dict = {}

    FORWARD_STRAND_FILTER = "fwd"
    REVERSE_STRAND_FILTER = "rev"
    NO_FILTER = "none"

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

    def set_geneset_filter(self, gene_set_file):
        # TODO the gene set file should be converted to a filter and then used to filter the transcript assembly
        # TODO (or should it be added to the mapping function?)
        # TODO (the time-efficient way of doing this is to generate a filtered list of transcripts and store them in memory: this list will be accessed many times afterwards)
        # TODO actually one can do that in the readers.py, when creating the annotation - only the transcripts from the given set will be retrieved from the transcript assembly and
        # TODO and stored in a list
        pass

    @preconditions(lambda span_size: isinstance(span_size, int))
    def yield_transcripts(self, span_size, filter=None):
        """
        Adds spanning regions to the start and the end of the transcript and yields the transcripts one-by-one.
        The original transcript assembly is left intact.

        :param filter: option to filter the transcripts
        :param span_size: int span_size - defines the number of nucleotides to span

        :return: iterable of type Transcript
        """

        if span_size < 0:
            error_message = "Negative span_size of %d provided: cannot span the transcripts with negative or " \
                            "zero values" % span_size
            self.logger.error(error_message)
            raise ValueError(error_message)
        self.logger.debug("Transcript span size: %d" % span_size)

        for transcript in self.transcript_assembly:
            span = transcript.spanning_segment

            new_start = span.start - span_size
            if new_start < 0:
                new_start = 0
            start_span = plastid.GenomicSegment(
                span.chrom, new_start, span.start, span.strand)

            # FIXME how do we check if we've reached the end?
            end_span = plastid.GenomicSegment(
                span.chrom, span.end, span.end + span_size, span.strand)

            spanned_transcript = transcript.__deepcopy__(None)
            spanned_transcript.add_segments(start_span, end_span)
            yield spanned_transcript

    # TODO in case gene_set filter is provided also filter with it

    @preconditions(lambda span_size: isinstance(span_size, int))
    def get_transcript_assembly(self, span_size, filter=None):
        """
        Returns transcripts if a transcript assembly already exists with provided span_size and filter.
        Creates and returns a new transcript assembly otherwise.

        :param span_size:  int span_size - defines the number of nucleotides to span
        :param gene_set_filter: list of genes to filter for
        :return: a vector of transcript or an iterable of type Transcript
        """

        if filter is None:
            filter = self.NO_FILTER

        if filter in self.transcript_assembly_dict and span_size in self.transcript_assembly_dict[filter]:
            return self.transcript_assembly_dict[filter][span_size]

        return self.generate_transcript_assembly(span_size, filter)

    @preconditions(lambda span_size: isinstance(span_size, int))
    def generate_transcript_assembly(self, span_size, filter=None):
        """
        Adds spanning regions to the start and the end of the transcript and yields the transcripts one-by-one.
        The original transcript assembly is left intact.

        :param span_size:  int span_size - defines the number of nucleotides to span
        :param gene_set_filter: list of genes to filter for

        :return: transcript vector
        """

        if span_size < 0:
            error_message = "Negative span_size of %d provided: cannot span the transcripts with negative or " \
                            "zero values" % span_size
            self.logger.error(error_message)
            raise ValueError(error_message)
        self.logger.debug("Transcript span size: %d" % span_size)

        if filter is not None:
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info("Generating transcript assembly under filter: %s and span size: %d"
                                                                 % (filter, span_size))
        this_transcript_assembly = [None]*self.transcript_count
        counter = 0
        for transcript in self.transcript_assembly:
            if filter == self.FORWARD_STRAND_FILTER:
                if transcript.strand != '+':
                    continue

            elif filter == self.REVERSE_STRAND_FILTER:
                if transcript.strand != '-':
                    continue

            span = transcript.spanning_segment

            # FIXME: I'm not sure this is a good solution: if you don't span where you have to,
            # FIXME: and the user assumes it is spanned this will lead to problems
            # FIXME: solution: return negative span size if allowed, and return 0's when retrieving counts on the negative segment

            new_start = span.start - span_size
            #if new_start < 0:
            #    new_start = 0
            start_span = plastid.GenomicSegment(
                span.chrom, new_start, span.start, span.strand)

            # FIXME how do we check if we've reached the end?
            end_span = plastid.GenomicSegment(
                span.chrom, span.end, span.end + span_size, span.strand)

            spanned_transcript = transcript.__deepcopy__(None)
            spanned_transcript.add_segments(start_span, end_span)
            this_transcript_assembly[counter] = spanned_transcript
            counter += 1

        if filter in self.transcript_assembly_dict:
            self.transcript_assembly_dict[filter].update({span_size:this_transcript_assembly})
        else:
            self.transcript_assembly_dict.update({filter:{span_size:this_transcript_assembly}})

        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info("Transcription assembly generation done")

        return this_transcript_assembly
    # TODO in case gene_set filter is provided also filter with it