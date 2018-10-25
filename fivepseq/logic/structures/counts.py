from math import floor

import plastid
from preconditions import preconditions
import numpy as np

from fivepseq import config
from logic.structures.alignment import Alignment
from logic.structures.annotation import Annotation


class FivePSeqCounts:
    START = "start"
    TERM = "termination"
    FULL_LENGTH = "full_length"

    alignment = None
    annotation = None

    counts_start = None
    counts_term = None
    counts_full_length = None
    meta_counts_start = None
    meta_counts_term = None

    def __init__(self, alignment, annotation):
        """
        Initializes a FivePSeqCounts object with Alignment and Annotation instances.
        :param alignment: Alignment type object
        :param annotation: Annotation type object
        """
        self.alignment = alignment
        self.annotation = annotation
        config.logger.debug("Initiated a FivePSeqCounts object with alignment from file %s and annotation from file %s."
                            % (alignment.alignment_file.filename, annotation.file_path))

    @preconditions(lambda span_size: isinstance(span_size, int),
                   lambda region: isinstance(region, str))
    def get_counts(self, span_size, region):
        """
        Returns arrays of read counts within regions spanning the given region, with the specified span_size.
        :param span_size: int: Specifies how many nucleotides to span around the specified region
        :param region: str: Specifies the region of the transcript to span for count vector generation
        :return: [[int]]: array of counts arrays of 5' mapping counts per position of the specified region of each transcript
        """

        if region == self.FULL_LENGTH:
            if self.counts_full_length is not None:
                return self.counts_full_length
        elif region == self.START:
            if self.counts_start is not None:
                return self.counts_start
        elif region == self.TERM:
            if self.counts_term is not None:
                return self.counts_term
        else:
            error_message = "Problem retrieving the counts. " \
                            "Invalid region \"%s\" specified: should be one of (%s, %s, %s)." \
                            % (region, self.FULL_LENGTH, self.START, self.TERM)
            config.logger.error(error_message)
            raise ValueError(error_message)

        config.logger.info("Retrieving counts from transcripts (span size :%d, position: %s)...."
                           % (span_size, region))

        count = self.annotation.transcript_count
        count_vectors = [None] * count
        progress_bar = ""
        i = 1
        transcript_generator = self.annotation.yield_transcripts(span_size)

        for transcript in transcript_generator:
            if i % 100 == 0:
                progress_bar += "#"
                print "\r>>Transcript count: %d (%d%s)\t%s" % (i, floor(100 * (i - 1) / count), '%', progress_bar),
            try:
                count_vector = self.get_transcript_counts(transcript, span_size, region)
            except Exception as e:
                raise e
            count_vectors[i - 1] = count_vector
            i += 1
        config.logger.info("Finished retrieving counts")
        return count_vectors

    @preconditions(lambda transcript: isinstance(transcript, plastid.genomics.roitools.Transcript),
                   lambda span_size: isinstance(span_size, int),
                   lambda region: isinstance(region, str))
    def get_transcript_counts(self, transcript, span_size, region):
        """
        Returns the array of counts per each transcript position within the given spanning region.

        :param region: str: Specifies the region of the transcript to span for count vector generation
        :param transcript: plastid.Transcript: The transcript to return the counts for: is is already spanned with the specified span_size
        :param span_size: int: Specifies how many nucleotides to span around the specified region
        :return: [int]: array of 5' mapping counts per position of the specified transcript region
        """

        try:
            count_vector = transcript.get_counts(self.alignment.bam_array)
            if region == self.FULL_LENGTH:
                pass
            elif region == self.START:
                count_vector = count_vector[:2 * span_size]
            elif region == self.TERM:
                count_vector = count_vector[-(2 * span_size):]
            else:
                error_message = "Problem retrieving the count vector for the transcript %s. " \
                                "Invalid region \"%s\" specified: should be one of (%s, %s, %s)." \
                                % (transcript.get_name(), region, self.FULL_LENGTH, self.START, self.TERM)
                config.logger.error(error_message)
                raise ValueError(error_message)
        except Exception as e:
            error_message = "Problem retrieving the count vector for the transcript %s. Reason:%s" % (
                transcript.get_name(), e.message)
            config.logger.error(error_message)
            raise Exception(error_message)

        count_vector = map(int, count_vector.tolist())

        return count_vector

    @preconditions(lambda span_size: isinstance(span_size, int),
                   lambda region: isinstance(region, str))
    def get_meta_counts(self, span_size, region):
        """
        Computes counts of 5' mapping positions at all the transcripts on the specified region, within the specified span size,
        and returns the position-wise sum of counts as a single [int] array.

        :param span_size: int: the number of bases to span around the specified region
        :param region: str: the region of transcript (start (START) or terminus (TERM)) to span around
        :return: [int] vector of position-wise sum of transcript-specific counts
        """

        if region == self.FULL_LENGTH:
            error_message = "Cannot compute meta counts for full length transcript counts: the counts should be of the same length. " \
                            "Regions can be specified from choices (%s, %s)" % (self.START, self.TERM)
            config.logger.error(error_message)
            raise ValueError(error_message)
        elif region == self.START:
            if self.meta_counts_start is not None:
                return self.meta_counts_start
        elif region == self.TERM:
            if self.meta_counts_term is not None:
                return self.meta_counts_term
        else:
            error_message = "Problem retrieving meta_counts. " \
                            "Invalid region \"%s\" specified: should be one of (%s, %s)." \
                            % (region, self.START, self.TERM)
            config.logger.error(error_message)
            raise ValueError(error_message)
        try:
            counts = self.get_counts(span_size, region)
        except Exception as e:
            raise e
        meta_counts = np.vstack(counts).sum(axis=0)
        if region == self.START:
            self.meta_counts_start = meta_counts
        elif region == self.TERM:
            self.meta_counts_term = meta_counts
        return meta_counts

    @staticmethod
    @preconditions(lambda counts: isinstance(counts, list),
                   lambda counts: isinstance(counts[0], int))
    def get_count_frames(counts):
        """
        Takes a vector of position-wise int counts and returns counts for three different frames.

        :param counts: [int]: a vector of int counts per-base
        :return: a tuple of frame count arrays (frame0:[int], frame1:[int], frame2:[int])
        """

        frame0_array = counts[range(0, len(counts), 3)]
        frame1_array = counts[range(1, len(counts), 3)]
        frame2_array = counts[range(2, len(counts), 3)]
        return frame0_array, frame1_array, frame2_array
