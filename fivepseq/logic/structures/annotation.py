"""
This module is meant to keep properties of annotation files and functions associated with those in one place .
"""
from preconditions import preconditions

import fivepseq.config
import plastid


class Annotation:
    file_path = None
    transcript_assembly = None
    geneset_filter = None
    transcript_count = None

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
            fivepseq.config.logger.error(error_message)
            raise ValueError(error_message)
        self.transcript_assembly = transcript_assembly
        self.transcript_count = len(transcript_assembly)
        fivepseq.config.logger.debug("Annotation object successfully created from file %s" % file_path)

    def set_geneset_filter(self, gene_set_file):
        # TODO the gene set file should be converted to a filter and then used to filter the transcript assembly
        # TODO (or should it be added to the mapping function?)
        # TODO (the time-efficient way of doing this is to generate a filtered list of transcripts and store them in memory: this list will be accessed many times afterwards)
        # TODO actually one can do that in the readers.py, when creating the annotation - only the transcripts from the given set will be retrieved from the transcript assembly and
        # TODO and stored in a list
        pass

    @preconditions(lambda span_size: isinstance(span_size, int))
    def yield_transcripts(self, span_size, gene_set_filter=None):
        """
        Adds spanning regions to the start and the end of the transcript and yields the transcripts one-by-one.
        The original transcript assembly is left intact.
        :param gene_set_filter: list of genes to filter for
        :param span_size: int span_size - defines the number of nucleotides to span
        :return: iterable of type Transcript
        """

        if span_size <= 0:
            error_message = "Negative or 0 span_size of %d provided: cannot span the transcripts with negative or " \
                            "zero values" % span_size
            fivepseq.config.logger.error(error_message)
            raise ValueError(error_message)
        fivepseq.config.logger.debug("Transcript span size: %d" % span_size)

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
