from math import floor

import plastid
import numpy as np
import pandas as pd

from preconditions import preconditions
from fivepseq import config


class FivePSeqCounts:
    """
    This class wraps annotation, alignment and genome objects in one place.
    Algorithms extracting count information from these objects are implemented in this class as functions.
    Algorithms able to work with count arrays and dataframes alone are in the algorithms package.
    """

    START = "start"
    TERM = "termination"
    FULL_LENGTH = "full_length"

    alignment = None
    annotation = None
    genome = None
    span_size = None

    count_vector_list_start = None
    count_vector_list_term = None
    count_vector_list_full_length = None
    meta_count_series_start = None
    meta_count_series_term = None

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

        config.logger.debug("Initiated a FivePSeqCounts object with"
                            "\n\talignment from file %s"
                            "\n\tannotation from file %s "
                            "\n\tgenome from file %s"
                            % (alignment.alignment_file.filename, annotation.file_path, genome.fasta_file))

    @preconditions(lambda region: isinstance(region, str))
    def get_count_vector_list(self, region):
        """
        Returns arrays of read count vectors spanning the given region of each transcript in the transcript assembly.
        The region is spanned according to the span_size set in self.

        :param region: str: Specifies the region of the transcript to span around
        :return: [[int]]: array of counts arrays of 5' mapping counts per position of the specified region of each transcript
        """

        # if counts are already computed, return the existing ones
        if region == self.FULL_LENGTH:
            if self.count_vector_list_full_length is not None:
                return self.count_vector_list_full_length
        elif region == self.START:
            if self.count_vector_list_start is not None:
                return self.count_vector_list_start
        elif region == self.TERM:
            if self.count_vector_list_term is not None:
                return self.count_vector_list_term
        else:
            error_message = "Cannot retrieve the counts. " \
                            "Invalid region \"%s\" specified: should be one of (%s, %s, %s)." \
                            % (region, self.FULL_LENGTH, self.START, self.TERM)
            config.logger.error(error_message)
            raise ValueError(error_message)

        # otherwise, retrieve the counts from the alignment file, referencing the transcript assembly
        config.logger.info("Retrieving counts (span size :%d, position: %s)..."
                           % (self.span_size, region))

        # initialize empty vectors
        count_vector_list = [None] * self.annotation.transcript_count

        # setup the progress bar and the counter
        progress_bar = ""
        counter = 1

        # loop through the transcripts yielded by the annotation object
        transcript_generator = self.annotation.yield_transcripts(self.span_size)
        for transcript in transcript_generator:

            # update progress to console
            if counter % 100 == 0:
                progress_bar += "#"
                config.logger.info("\r>>Transcript count: %d (%d%s)\t%s" % (
                    counter, floor(100 * (counter - 1) / self.annotation.transcript_count), '%', progress_bar), )

            # retrieve actual counts for current transcript
            try:
                count_vector = self.get_count_vector(transcript, self.span_size, region)
                count_vector_list[counter - 1] = count_vector

            except Exception as e:
                error_message = "Problem retrieving counts for transcript %s. Reason: %s" \
                                % (transcript.get_name, e.message)
                config.logger.error(error_message)
                raise Exception(error_message)

            counter += 1

        # save the list for further reference
        self.set_count_vector_list(count_vector_list, region)

        # report successful retrieval
        config.logger.info("Finished retrieving count vectors")

        return count_vector_list

    @preconditions(lambda transcript: isinstance(transcript, plastid.genomics.roitools.Transcript),
                   lambda span_size: isinstance(span_size, int),
                   lambda region: isinstance(region, str))
    def get_count_vector(self, transcript, span_size, region):
        """
        Returns the vector of counts for the given transcript within the given spanning region.

        :param region: str: Specifies the region of the transcript to span for count vector generation
        :param transcript: plastid.Transcript: The transcript to return the counts for: is is already spanned with the specified span_size
        :param span_size: int: Specifies how many nucleotides to span around the specified region
        :return: [int]: array of 5' mapping counts per position of the specified transcript region
        """

        try:
            # retrieve the count vector using plastid function "get_counts" called from the given Transcript object
            count_vector = transcript.get_counts(self.alignment.bam_array)

            # return only the region of the vector that is specified by region and span_size parameters
            if region == self.FULL_LENGTH:
                # the full vector will be returned
                pass
            elif region == self.START:
                count_vector = count_vector[:2 * span_size]
            elif region == self.TERM:
                count_vector = count_vector[-(2 * span_size):]
            else:
                error_message = "Cannot retrieve a count vector for the transcript %s. " \
                                "Invalid region \"%s\" specified: should be one of (%s, %s, %s)." \
                                % (transcript.get_name(), region, self.FULL_LENGTH, self.START, self.TERM)
                config.logger.error(error_message)
                raise ValueError(error_message)

        except Exception as e:
            error_message = "Problem retrieving the count vector for the transcript %s. Reason:%s" % (
                transcript.get_name(), e.message)
            config.logger.error(error_message)
            raise Exception(error_message)

        # convert the count array to an int vector
        count_vector = map(int, count_vector.tolist())

        return count_vector

    @preconditions(lambda region: isinstance(region, str))
    def get_meta_count_series(self, region):
        """
        Computes counts of 5' mapping positions at all the transcripts on the specified region, within the specified span size,
        and returns the position-wise sum of counts as a single [int] array.

        :param region: str: the region of transcript (start (START) or terminus (TERM)) to span around

        :return: pd.Series{int: int}: series of position-wise sum of transcript-specific counts indexed according to the
        distance of genomic coordinates from the first nucleotides of the codon corresponding to the specified region (START or TERM)
        """

        if region == self.FULL_LENGTH:
            error_message = "Cannot compute meta counts for full length transcript counts: the counts should be of " \
                            "the same length. " \
                            "Regions can be specified from choices (%s, %s)" % (self.START, self.TERM)
            config.logger.error(error_message)
            raise ValueError(error_message)
        elif region == self.START:
            if self.meta_count_series_start is not None:
                return self.meta_count_series_start
        elif region == self.TERM:
            if self.meta_count_series_term is not None:
                return self.meta_count_series_term
        else:
            error_message = "Problem retrieving meta_counts. " \
                            "Invalid region \"%s\" specified: should be one of (%s, %s)." \
                            % (region, self.START, self.TERM)
            config.logger.error(error_message)
            raise ValueError(error_message)
        try:
            count_vector_list = self.get_count_vector_list(region)
        except Exception as e:
            raise e
        meta_count_series = CountManager.meta_count_vector_to_series(
            CountManager.compute_meta_counts(count_vector_list), region, tail=self.span_size)

        self.set_meta_count_series(meta_count_series, region)

        return meta_count_series

    @preconditions(lambda count_vector_list: isinstance(count_vector_list, list),
                   lambda count_vector_list: isinstance(count_vector_list[0], list),
                   lambda count_vector_list: isinstance(count_vector_list[0][0], int),
                   lambda region: isinstance(region, str))
    def set_count_vector_list(self, count_vector_list, region):
        """
        Sets the retrieved counts as a class property for later use. The property is chosen is based on the region supplied.

        :param count_vector_list: [[int]]: the vector of count vectors per transcript
        :param region: str: the region for which the counts were computed

        :return: nothing to return
        """
        if region == self.START:
            self.count_vector_list_start = count_vector_list
        elif region == self.TERM:
            self.count_vector_list_term = count_vector_list
        elif region == self.FULL_LENGTH:
            self.count_vector_list_full_length = count_vector_list
        else:
            error_message = "Cannot set counts: wrong region %s supplied: should be either of (%s, %s, %s)" \
                            % (region, self.START, self.TERM, self.FULL_LENGTH)
            config.logger.error(error_message)
            raise ValueError(error_message)

    @preconditions(lambda meta_count_series: isinstance(meta_count_series, pd.Series),
                   lambda region: isinstance(region, str))
    def set_meta_count_series(self, meta_count_series, region):
        """
        Sets the retrieved meta-counts as a class property for later use. The property is chosen is based on the region supplied.

        :param meta_count_series: Series{int:int}: the panda Series of per-position mapped read sums across transcripts
        indexed by position from first nucleotide of START of STOP codon
        :param region: str: the region for which the counts were computed

        :return: nothing to return
        """
        if region == self.START:
            self.meta_count_series_start = meta_count_series
        elif region == self.TERM:
            self.meta_count_series_term = meta_count_series

    @preconditions(lambda region: isinstance(region, str),
                   lambda span_before: isinstance(span_before, int),
                   lambda span_before: span_before >= 0,
                   lambda span_after: isinstance(span_after, int),
                   lambda span_after: span_after >= 0)
    def get_unique_sequences(self, region, span_before, span_after):
        """
        Retrieves the unique set of sequences spanning the given region of all transcripts, with the specified parameters.

        :param region: str: the START or TERM parts of the transcript
        :param span_before: int: the number of nucleotides to span before the first codon of the specified region
        :param span_after: int: the number of nucleotides to span after the last codon of the specified region

        :return: dict{str:int}: a dictionary keyed by the unique sequences identified within the spanning regions
        and valued by the number of occurrences of that sequence in transcripts
        """
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


class CountManager:
    """
    This module impelements a set of static functions to handle count vectors retrieved from FivePSeqCounts class.
    """

    def __init__(self):
        pass

    @staticmethod
    @preconditions(lambda count_vector_list: isinstance(count_vector_list, list),
                   lambda count_vector_list: isinstance(count_vector_list[0], list),
                   lambda count_vector_list: isinstance(count_vector_list[0][0], int))
    def compute_meta_counts(count_vector_list):
        """
        Computes the sum of counts at each position across transcripts.

        :param count_vector_list: [[int]] a list of count vectors for all transcripts

        :return: [int]: a vector of position-wise count sums
        """

        # TODO check that the count vectors have the same length

        # sum the position-wise counts
        meta_count_vector = np.vstack(count_vector_list).sum(axis=0).tolist()

        return meta_count_vector

    @staticmethod
    @preconditions(lambda count_vector: isinstance(count_vector, list),
                   lambda count_vector: isinstance(count_vector[0], int),
                   lambda span_size: isinstance(span_size, int),
                   lambda region: isinstance(region, str),
                   lambda include_span: isinstance(include_span, bool))
    def extract_frame_count_vectors(count_vector, span_size, region=FivePSeqCounts.START, include_span=False):
        """
        Takes a vector of position-wise int counts and returns counts for three different frames from 0 to 2,
        relative to either START (default) or TERM regions.

        :param count_vector: [int]: a transcript-specific vector of 5P read counts per transcript position
        :param span_size: int: the size of regions spanning around the transcript cds
        :param include_span: if true returns frame counts including the spanning regions, and only the cds region
        otherwise
        :param region: region (START or TERM) relative to which to count the frames

        :return: a tuple of frame count arrays (frame0:[int], frame1:[int], frame2:[int])
        """

        # determine the tail size to be subtracted from the count_vector
        if include_span:
            tail = 0
        else:
            tail = span_size

        # for START, start the Frame0 from tail to the length of the vector minus the tail
        if region == FivePSeqCounts.START:
            frame0_array = count_vector[0 + tail: len(count_vector) - tail: 3]
            frame1_array = count_vector[1 + tail: len(count_vector) - tail: 3]
            frame2_array = count_vector[2 + tail: len(count_vector) - tail: 3]

        # for TERM, align the Frame0 to the last nucleotide of the cds (exclude span, tail = span_size)
        # or the count_vector (include span, tail = 0)
        elif region == FivePSeqCounts.TERM:
            # NOTE the frames relative to START and TERM should be aligned in the future
            # NOTE (if cds length is not a multiple of 3)
            frame0_array = [count_vector[i] for i in list(reversed(range(len(count_vector) - 1 - tail, -1 + tail, -3)))]
            frame1_array = [count_vector[i] for i in list(reversed(range(len(count_vector) - 2 - tail, -1 + tail, -3)))]
            frame2_array = [count_vector[i] for i in list(reversed(range(len(count_vector) - 3 - tail, -1 + tail, -3)))]

        else:
            error_message = "Invalid region %s specified: should be either %s or %s" \
                            % (region, FivePSeqCounts.START, FivePSeqCounts.TERM)
            config.logger.error(error_message)
            raise Exception(error_message)

        return frame0_array, frame1_array, frame2_array

    @staticmethod
    @preconditions(lambda meta_count_vector: isinstance(meta_count_vector, list),
                   lambda meta_count_vector: isinstance(meta_count_vector[0], int),
                   lambda region: isinstance(region, str),
                   lambda tail: isinstance(tail, int),
                   lambda tail: tail >= 0)
    def meta_count_vector_to_series(meta_count_vector, region, tail=0):
        """
        Takes a vector of meta counts, indexes them with distances from the specified region.
        Returns a series with indexes as genomic coordinates from start/stop codons, and values as meta counts at each coordinates.

        For the following REAL coordinates (R), D0, D1 and D2 will be converted to:

            Relative to START:

            R:  0,      1,                  ... S0,     Sm,     S1,     ... L - 1
            D: -T,     -T+1,                ... 0,      1,      2,      ... L - T - 1

            Relative to TERM:

            R:  0,          1,              ... S0,     Sm,     S1,     ... L - 1
            D: -(L-T-3),     -(L-T-2),      ... 0,      1,      2,      ... T + 2


            Legend:
                R   Real coordinates
                S0  First nucleotides of START or TERM codon
                S1  Last nucleotides of START or TERM codon
                T   Tail
                L   Vector length

        :param meta_count_vector: [int]: a vector of summed position-wise counts for a meta-transcript
        :param region: str: the region respective to which the distance is calculated
        :param tail: int:

        :return: pandas.Series: a series with indices as genomic coordinates* and values as meta counts.
        *-corresponding to positions' distance from nucleotide 0 of START/TERM codons

        """

        if region == FivePSeqCounts.START:
            d = np.arange(-tail, len(meta_count_vector) - tail )

        elif region == FivePSeqCounts.TERM:
            d = np.arange(-(len(meta_count_vector) - tail - 3), tail + 3)

        else:
            error_message = "Invalid region %s specified: should be either %s or %s" \
                            % (region, FivePSeqCounts.START, FivePSeqCounts.TERM)
            config.logger.error(error_message)
            raise Exception(error_message)

        meta_counts_series = pd.Series(data=meta_count_vector, index=d)

        return meta_counts_series

    @staticmethod
    @preconditions(lambda region: isinstance(region, str),
                   lambda span_size: isinstance(span_size, int))
    def extract_count_sums_per_frame_per_transcript(count_vector_list, span_size, region):
        """
        Returns a data frame with rows representing transcripts and columns (F0, F1, F2) representing the sum of 5P
        read mapping counts at each frame. The transcripts are aligned at the start or the end, depending on the
        region specified.

        :param span_size:
        :param count_vector_list: the list of per-transcript count vectors
        :param region: str: the region to align the transcripts to

        :return: a dataframe with frame-based count-sums for each transcript
        """

        config.logger.debug("Retrieving count-sums per frame relative to %s ..." % region)

        # Create an empty dataframe
        n = len(count_vector_list)
        frame_counts_df = pd.DataFrame({'F0': [0] * n,
                                        'F1': [0] * n,
                                        'F2': [0] * n})
        progress_bar = ""
        for t_ind in range(0, n):
            # Print status update to console
            if t_ind % 100 == 0:
                progress_bar += "#"
                print "\r>>Transcript count: %d (%d%s)\t%s" % (
                    t_ind, np.floor(100 * (t_ind - 1) / n), '%', progress_bar),

            # extract frame count vectors from count vectors
            count_vector = count_vector_list[t_ind]
            frame_counts = CountManager.extract_frame_count_vectors(count_vector, span_size, region)

            # sum-up counts in each frame and add to the dataframe
            for f in range(0, 3):
                frame_counts_df.iloc[t_ind, f] = sum(frame_counts[f])

        return frame_counts_df
