from math import floor

import plastid
from preconditions import preconditions
import numpy as np
import pandas as pd

from fivepseq import config
from logic.structures.alignment import Alignment
from logic.structures.annotation import Annotation


class FivePSeqCounts:
    START = "start"
    TERM = "termination"
    FULL_LENGTH = "full_length"

    alignment = None
    annotation = None
    genome = None
    span_size = None

    counts_start = None
    counts_term = None
    counts_full_length = None
    meta_counts_start = None
    meta_counts_term = None

    check_stop_codons = True

    @preconditions(lambda span_size: isinstance(span_size, int))
    def __init__(self, alignment, annotation, genome, span_size):
        """
        Initializes a FivePSeqCounts object with Alignment and Annotation instances.
        :param alignment: fivepseq.logic.structures.Alignment type object
        :param annotation: fivepseq.logic.structures.Annotation type object
        :param genome: fivepseq.logic.structures.Genome: Genome type object
        :param span_size: int: Specifies how many nucleotides to span around any specified region of transcripts
        """
        self.alignment = alignment
        self.annotation = annotation
        self.genome = genome
        self.span_size = span_size
        config.logger.debug("Initiated a FivePSeqCounts object with alignment from file %s and annotation from file %s."
                            % (alignment.alignment_file.filename, annotation.file_path))

    @preconditions(lambda region: isinstance(region, str))
    def get_counts(self, region):
        """
        Returns arrays of read counts within regions spanning the given region, with the specified span_size.
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
                           % (self.span_size, region))

        count = self.annotation.transcript_count
        count_vectors = [None] * count
        progress_bar = ""
        i = 1
        transcript_generator = self.annotation.yield_transcripts(self.span_size)

        for transcript in transcript_generator:
            if i % 100 == 0:
                progress_bar += "#"
                print "\r>>Transcript count: %d (%d%s)\t%s" % (i, floor(100 * (i - 1) / count), '%', progress_bar),
            try:
                count_vector = self.get_transcript_counts(transcript, self.span_size, region)
            except Exception as e:
                raise e
            count_vectors[i - 1] = count_vector
            i += 1
        config.logger.info("Finished retrieving counts")

        self.set_counts(count_vectors, region)
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

    @preconditions(lambda region: isinstance(region, str))
    def get_meta_counts(self, region):
        """
        Computes counts of 5' mapping positions at all the transcripts on the specified region, within the specified span size,
        and returns the position-wise sum of counts as a single [int] array.

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
            counts = self.get_counts(region)
        except Exception as e:
            raise e
        meta_counts = np.vstack(counts).sum(axis=0).tolist()

        self.set_meta_counts(meta_counts, region)

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

        frame0_array = counts[0: len(counts): 3]
        frame1_array = counts[1: len(counts): 3]
        frame2_array = counts[2: len(counts): 3]
        return frame0_array, frame1_array, frame2_array

    @preconditions(lambda region: isinstance(region, str))
    def get_meta_counts_df(self, region):
        #TODOC
        try:
            meta_counts = self.get_meta_counts(region)
        except Exception as e:
            raise e

        distance_from_last = np.arange(-self.span_size, self.span_size)
        distance_from_first = distance_from_last + 3

        df = pd.DataFrame({'D1': distance_from_first,
                           'D2': distance_from_last,
                           'Count': meta_counts})
        return df

    @preconditions(lambda count_vectors: isinstance(count_vectors, list),
                   lambda count_vectors: isinstance(count_vectors[0], list),
                   lambda count_vectors: isinstance(count_vectors[0][0], int),
                   lambda region: isinstance(region, str))
    def set_counts(self, count_vectors, region):
        """
        Sets the retrieved counts as a class property for later use. The property is chosen is based on the region supplied.
        :param count_vectors: [[int]]: the vector of count vectors per transcript
        :param region: str: the region for which the counts were computed
        :return: nothing to return
        """
        if region == self.START:
            self.counts_start = count_vectors
        elif region == self.TERM:
            self.counts_term = count_vectors
        elif region == self.FULL_LENGTH:
            self.counts_full_length = count_vectors
        else:
            error_message = "Cannot set counts: wrong region %s supplied: should be either of (%s, %s, %s)" \
                            % (region, self.START, self.TERM, self.FULL_LENGTH)
            config.logger.error(error_message)
            raise ValueError(error_message)

    # @preconditions(lambda meta_counts: isinstance(meta_counts, list),
    #               lambda meta_counts: isinstance(meta_counts[0], int),
    #               lambda region: isinstance(region, str))
    def set_meta_counts(self, meta_counts, region):
        """
        Sets the retrieved meta-counts as a class property for later use. The property is chosen is based on the region supplied.
        :param meta_counts: [int]: the vector of per-position mapped read sums across transcripts
        :param region: str: the region for which the counts were computed
        :return: nothing to return
        """
        if region == self.START:
            self.meta_counts_start = meta_counts
        elif region == self.TERM:
            self.meta_counts_term = meta_counts

    def get_unique_sequences(self, region, span_before, span_after):
        sequences = {}
        i = 0
        for transcript in self.annotation.yield_transcripts(max(span_before, span_after)):
            sequence = transcript.get_sequence(self.genome.genome_dict)
            if region == self.TERM:
                endpoint = len(transcript.spanning_segment) - max(span_before, span_after)
                span_sequence = sequence[endpoint - span_before: endpoint + span_after]
            elif region == self.START:
                startpoint = max(span_before, span_after)
                span_sequence = sequence[startpoint - span_before: startpoint + span_after]
            else:
                raise Exception
            if span_sequence in sequences.keys():
                sequences[span_sequence] += 1
            else:
                sequences[span_sequence] = 1
            i += 1
        return sequences
