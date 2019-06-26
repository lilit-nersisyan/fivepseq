import logging
import os
from math import floor

import numpy as np
import pandas as pd
import plastid
from preconditions import preconditions
from scipy import stats

from fivepseq import config
from fivepseq.logic.structures import codons
from fivepseq.logic.structures.codons import Codons
from fivepseq.util.writers import FivePSeqOut


class FivePSeqCounts:
    """
    This class wraps annotation, alignment and genome objects in one place.
    Algorithms extracting count information from these objects are implemented in this class as functions.
    Algorithms able to work with count arrays and dataframes alone are in the algorithms package.
    """

    START = "start"
    TERM = "termination"
    FULL_LENGTH = "full_length"
    ALL = "all"

    START_CODON = "start"
    STOP_CODON = "stop"
    TRANSCRIPT_LENGTH = "len"
    TRANSCRIPT_3NT = "3nt"
    NUMBER_READS = "NumOfReads"
    NUMBER_POSITIONS = "NumOfMapPositions"

    count_distribution = None
    outlier_lower = None
    downsample_by = None

    alignment = None
    annotation = None
    genome = None

    count_vector_list_start = None
    count_vector_list_term = None
    count_vector_list_full_length = None
    meta_count_series_start = None
    meta_count_series_term = None
    frame_counts_df_start = None
    frame_counts_df_term = None

    outliers = []

    start_codon_dict = {}
    stop_codon_dict = {}
    canonical_transcript_index = []
    transcript_descriptors = None

    loci_file = None

    logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)

    def __init__(self, alignment, annotation, genome, downsample_constant, transcript_filter=None):
        """
        Initializes a FivePSeqCounts object with Alignment and Annotation instances.

        :param alignment: fivepseq.logic.structures.Alignment type object
        :param annotation: fivepseq.logic.structures.Annotation type object
        :param genome: fivepseq.logic.structures.Genome: Genome type object
        """

        self.alignment = alignment
        self.annotation = annotation
        self.genome = genome
        self.transcript_filter = transcript_filter
        self.downsample_by = downsample_constant
        self.logger.info("Initiated a FivePSeqCounts object with"
                         "\n\talignment from file %s"
                         "\n\tannotation from file %s "
                         "\n\tgenome from file %s"
                         % (alignment.alignment_file.filename, annotation.file_path, genome.fasta_file))

    def get_transcript_descriptors(self):
        if self.transcript_descriptors is None:
            self.generate_transcript_descriptors()

        return self.transcript_descriptors

    def get_start_codon_dict(self):
        if self.start_codon_dict is None:
            self.generate_transcript_descriptors()

        return self.start_codon_dict

    def get_stop_codon_dict(self):
        if self.stop_codon_dict is None:
            self.generate_transcript_descriptors()

        return self.stop_codon_dict

    def generate_transcript_descriptors(self):
        """
        Generates and stores the basic statistics on transcript sequences and counts.
        The following objects are generated and kept in self:

        transcript_descriptors:: pandas DataFrame
            - columns: START, TERM codons, transcript length,
                transcript length divisible by three, number of reads mapping within coding region
            - rows: transcripts


        :return:
        """

        # info
        self.logger.info("Generating transcript descriptors")

        transcript_assembly = self.annotation.get_transcript_assembly_default_filter(0)
        transcript_count = len(transcript_assembly)
        self.transcript_descriptors = pd.DataFrame(data=None,
                                                   index=range(transcript_count),
                                                   columns=[self.START_CODON,
                                                            self.STOP_CODON,
                                                            self.TRANSCRIPT_LENGTH,
                                                            self.TRANSCRIPT_3NT,
                                                            self.NUMBER_READS,
                                                            self.NUMBER_POSITIONS])

        self.count_distribution = []

        for transcript_ind in range(transcript_count):
            transcript = transcript_assembly[transcript_ind]

            cds_sequence = self.get_cds_sequence_safe(transcript, 0)
            count_vector = self.get_count_vector_safe(transcript, 0)

            # NOTE the count distribution does not include values 0 and 1 to avoid skewness for outlier detection
            for c in count_vector:
                if c > 0:
                    self.count_distribution.append(c)

            start_codon = cds_sequence[0:3]
            stop_codon = cds_sequence[len(cds_sequence) - 3:len(cds_sequence)]

            if (start_codon == codons.Codons.START_CODON) & (stop_codon in codons.Codons.stop_codons):
                self.canonical_transcript_index.append(transcript_ind)

            self.transcript_descriptors.at[transcript_ind, self.START_CODON] = start_codon
            self.transcript_descriptors.at[transcript_ind, self.STOP_CODON] = stop_codon
            self.transcript_descriptors.at[transcript_ind, self.TRANSCRIPT_3NT] = str(len(cds_sequence) % 3 == 0)
            self.transcript_descriptors.at[transcript_ind, self.TRANSCRIPT_LENGTH] = len(cds_sequence)
            self.transcript_descriptors.at[transcript_ind, self.NUMBER_READS] = int(np.sum(count_vector))
            self.transcript_descriptors.loc[transcript_ind, self.NUMBER_POSITIONS] = np.count_nonzero(count_vector)

            if start_codon in self.start_codon_dict.keys():
                self.start_codon_dict[start_codon] += 1
            else:
                self.start_codon_dict.update({start_codon: 1})

            if stop_codon in self.stop_codon_dict.keys():
                self.stop_codon_dict[stop_codon] += 1
            else:
                self.stop_codon_dict.update({stop_codon: 1})

        self.outlier_lower = self.get_outlier_lower()

        self.logger.info("The lower bound for outliers set as %f " % self.outlier_lower)

        self.logger.info("Done generating transcript descriptors")

    def get_count_distribution(self):
        if self.count_distribution is not None:
            return self.count_distribution

        else:
            self.generate_transcript_descriptors()
            return self.count_distribution

    @preconditions(lambda count_distribution: isinstance(count_distribution, list))
    def set_count_distribution(self, count_distribution):
        """
        Sets the count distribution according to the specified count vector.

        :param count_distribution: a vector of counts (should be [int] but [float] is also acceptable)
        :return:
        :raise: ValueError if values in the count distribution not convertable to int
        """
        if len(count_distribution) == 0:
            self.count_distribution = []
        else:
            try:
                count_distribution = map(int, count_distribution)
                self.count_distribution = count_distribution
            except Exception as e:
                raise ValueError("problem converting count distribution values to int: %s" % str(e))

    def get_outlier_lower(self):
        """
        Returns the lower bound for outliers detected as points lying self.downsample_by number times higher than the
        25-75% interquartile range.

        :param count_distribution:
        :return:
        """
        if self.outlier_lower is not None:
            return self.outlier_lower

        count_distribution = self.get_count_distribution()
        if len(count_distribution) == 0:
            self.outlier_lower = 0
            return 0

        scd = sorted(count_distribution)
        lam = np.mean(scd)
        ps = [1 - stats.poisson.cdf(x, lam) for x in scd]
        ind = np.where(np.asarray(ps) <= 0)[0].tolist()

        if len(ind) > 0:
            # outliers = [scd[i] for i in ind]
            outlier_lower = scd[min(ind) - 1]
        else:
            outlier_lower = max(scd) + 1

        self.outlier_lower = outlier_lower
        return outlier_lower

    def set_outlier_lower(self, outlier_lower):
        # TODO add checks, set preconditions
        self.outlier_lower = outlier_lower

    def generate_count_vector_lists(self):
        """
        Generates read count vectors for full length transcripts, terminus- and start- aligned sections,
        spanning respective regions of each transcript in the transcript assembly.
        The region is spanned according to the span_size set in annotation.

        :return: [[int]]: array of counts arrays of 5' mapping counts per position of the specified region of each transcript
        """

        # if counts are already computed, return the existing ones

        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info("Generating count vectors")

        if self.count_vector_list_full_length is not None:
            if self.count_vector_list_term is not None:
                if self.count_vector_list_start is not None:
                    logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).warning("All count vectors are already generated")

        # otherwise, retrieve the counts from the alignment file, referencing the transcript assembly
        self.logger.info("Retrieving counts (span size :%d)..."
                         % self.annotation.span_size)

        # initialize empty vectors
        transcript_count = len(self.annotation.get_default_transcript_assembly())
        self.count_vector_list_full_length = [None] * transcript_count
        self.count_vector_list_term = [None] * transcript_count
        self.count_vector_list_start = [None] * transcript_count

        # setup the the counter
        counter = 1

        for transcript in self.annotation.get_default_transcript_assembly():

            if transcript.get_name() == "rna-XM_021479082.1":
                print "here"

            # update to console
            if counter % 10000 == 0:
                self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                    counter, floor(100 * (counter - 1) / self.annotation.transcript_count),
                    '%'), )

            # retrieve actual counts for current transcript
            try:
                count_vector = self.get_count_vector(transcript, self.annotation.span_size, self.FULL_LENGTH)
                self.count_vector_list_full_length[counter - 1] = count_vector
                self.count_vector_list_start[counter - 1] = count_vector[:2 * self.annotation.span_size]
                self.count_vector_list_term[counter - 1] = count_vector[-(2 * self.annotation.span_size):]

            except Exception as e:
                error_message = "Problem retrieving counts for transcript %s. Reason: %s" \
                                % (transcript.get_name(), e.message)
                self.logger.error(error_message)
                raise Exception(error_message)

            counter += 1
        self.check_for_codons = False

        # report successful retrieval
        self.logger.info("Finished retrieving count vectors")

    @preconditions(lambda region: isinstance(region, str))
    def get_count_vector_list(self, region):
        """
        Returns arrays of read count vectors spanning the given region of each transcript in the transcript assembly.
        The region is spanned according to the span_size set in annotation.

        :param region: str: Specifies the region of the transcript to span around
        :return: [[int]]: array of counts arrays of 5' mapping counts per position of the specified region of each transcript
        """

        # if counts are already computed, return the existing ones else generate count vector lists first
        if self.count_vector_list_full_length is None:
            self.generate_count_vector_lists()

        if region == self.FULL_LENGTH:
            return self.count_vector_list_full_length
        elif region == self.START:
            return self.count_vector_list_start
        elif region == self.TERM:
            return self.count_vector_list_term

        else:
            error_message = "Cannot retrieve the counts. " \
                            "Invalid region \"%s\" specified: should be one of (%s, %s, %s)." \
                            % (region, self.FULL_LENGTH, self.START, self.TERM)
            self.logger.error(error_message)
            raise ValueError(error_message)

    @preconditions(lambda transcript: isinstance(transcript, plastid.genomics.roitools.Transcript),
                   lambda span_size: isinstance(span_size, int),
                   lambda region: isinstance(region, str))
    def get_count_vector(self, transcript, span_size, region, downsample=True):
        """
        Returns the vector of counts for the given transcript within the given spanning region.

        :param region: str: Specifies the region of the transcript to span for count vector generation
        :param transcript: plastid.Transcript: The transcript to return the counts for: is is already spanned with the specified span_size
        :param span_size: int: Specifies how many nucleotides to span around the specified region
        :param transcript_ind: int: the index of transcript in the transcript assembly
        :return: [int]: array of 5' mapping counts per position of the specified transcript region
        """

        try:
            # retrieve the count vector using plastid function "get_counts" called from the given Transcript object
            count_vector = self.get_count_vector_safe(transcript, span_size)
            if downsample and np.array(count_vector > self.outlier_lower).any():
                count_vector_ds = [0] * len(count_vector)
                for i in range(len(count_vector_ds)):
                    if count_vector[i] > self.outlier_lower:
                        count_vector_ds[i] = self.outlier_lower
                        outlier_params = [FivePSeqOut.get_transcript_attr(transcript, "ID"),
                                          FivePSeqOut.get_transcript_attr(transcript, "Name"),
                                          i - span_size, len(count_vector) - i - span_size, count_vector[i],
                                          count_vector_ds[i]]
                        if outlier_params not in self.outliers:
                            self.outliers.append(outlier_params)
                    else:
                        count_vector_ds[i] = count_vector[i]
                count_vector = count_vector_ds

            count_vector = count_vector[transcript.cds_start: transcript.cds_end + 2 * span_size]

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
                self.logger.error(error_message)
                raise ValueError(error_message)

        except Exception as e:
            error_message = "Problem retrieving the count vector for the transcript %s. Reason:%s" % (
                transcript.get_name(), e.message)
            self.logger.error(error_message)
            raise Exception(error_message)

        # convert the count array to an int vector
        if not isinstance(count_vector, list):
            count_vector = count_vector.tolist()
        # if not isinstance(count_vector[0], int):
        count_vector = map(int, count_vector)

        return count_vector

    def get_count_vector_safe(self, transcript, span_size):
        """
        A safe method to return count vector accounting for transcripts that span before or after genome start and end.

        :param transcript:
        :param span_size:
        :return:
        """

        try:
            count_vector = transcript.get_counts(self.alignment.bam_array)
        except Exception as e:

            if transcript.spanning_segment.start < 0:
                diff = -1 * transcript.spanning_segment.start
                t_subchain = transcript.get_subchain(diff, transcript.spanning_segment.end, stranded=False)
                subchain_counts = list(t_subchain.get_counts(self.alignment.bam_array))
                count_vector = [0] * diff + subchain_counts

                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER). \
                    debug("Transcript %s at the beginning of the genome padded with %d zeros"
                          % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))

            else:
                t_len = transcript.spanning_segment.end - transcript.spanning_segment.start
                diff = transcript.spanning_segment.end - len(self.genome.genome_dict[transcript.chrom].seq)

                if diff > span_size:
                    # NOTE wrongly annotated transcripts go outside genome boundaries,
                    # NOTE return an empty vector spanned by span size as a safe way of discarding such transcripts
                    count_vector = [0] * t_len
                    logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER). \
                        debug("Transcript %s exceeds genome dimensions by %d bases"
                              % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))

                else:
                    t_subchain = transcript.get_subchain(diff, transcript.spanning_segment.end, stranded=False)

                    subchain_counts = list(t_subchain.get_counts(self.alignment.bam_array))
                    count_vector = subchain_counts + [0] * diff

                    logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER). \
                        debug("Transcript %s at the end of the genome padded with %d zeros"
                              % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))

        return count_vector

    def get_cds_sequence_safe(self, transcript, span_size):
        try:
            sequence = transcript.get_sequence(self.genome.genome_dict)
            cds_sequence = sequence[transcript.cds_start + span_size: transcript.cds_end + span_size]
        except:

            if transcript.spanning_segment.start < 0:
                diff = -1 * transcript.spanning_segment.start
                t_subchain = transcript.get_subchain(diff, transcript.spanning_segment.end, stranded=False)
                sequence = t_subchain.get_sequence(self.genome.genome_dict)

                if span_size < diff:
                    cds_sequence = sequence[transcript.cds_start + span_size - diff: transcript.cds_end + span_size]
                else:  # TODO I don't know how to get sequence in this case: need debugging
                    cds_sequence = sequence[transcript.cds_start + span_size - diff: transcript.cds_end + span_size]

                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER). \
                    debug("Transcript %s at the beginning of the genome padded with %d N's"
                          % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))

            else:
                t_len = transcript.spanning_segment.end - transcript.spanning_segment.start
                diff = transcript.spanning_segment.end - len(self.genome.genome_dict[transcript.chrom].seq)

                if diff > span_size:
                    # NOTE wrongly annotated transcripts go outside genome boundaries,
                    # NOTE return an empty sequence spanned by span size as a safe way of discarding such transcripts
                    cds_sequence = ''.join(['N'] * t_len)

                else:
                    t_subchain = transcript.get_subchain(diff, transcript.spanning_segment.end, stranded=False)

                    sequence = t_subchain.get_sequence(self.genome.genome_dict)
                    cds_sequence = sequence[transcript.cds_start + span_size:
                                            transcript.cds_end + span_size - diff]

        return cds_sequence

    def get_outliers_df(self):
        """
        Returns the outliers in the form of a data-frame with column names.

        :return:
        """

        colnames = ["ID", "Name", "position_from_start", "position_from_term", "actual_count", "downsampled_count"]
        outliers_df = pd.DataFrame(self.outliers, index=None, columns=colnames)
        return outliers_df

    @preconditions(lambda region: isinstance(region, str))
    def get_frame_counts_df(self, region):
        if region == self.START:
            if self.frame_counts_df_start is None:
                self.frame_counts_df_start = CountManager.extract_count_sums_per_frame_per_transcript(
                    self.get_count_vector_list(FivePSeqCounts.FULL_LENGTH), self.annotation.span_size,
                    FivePSeqCounts.START)
            return self.frame_counts_df_start

        elif region == self.TERM:
            if self.frame_counts_df_term is None:
                self.frame_counts_df_term = CountManager.extract_count_sums_per_frame_per_transcript(
                    self.get_count_vector_list(FivePSeqCounts.FULL_LENGTH), self.annotation.span_size,
                    FivePSeqCounts.TERM)
            return self.frame_counts_df_term

        else:
            err_msg = ("Wrong region %s provided: should be either %s or %s"
                       % (region, self.START, self.TERM))
            self.logger.error(err_msg)
            raise Exception(err_msg)

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
            self.logger.error(error_message)
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
            self.logger.error(error_message)
            raise ValueError(error_message)
        try:
            count_vector_list = self.get_count_vector_list(region)
        except Exception as e:
            raise e
        meta_count_series = CountManager.count_vector_to_series(
            CountManager.compute_meta_counts(count_vector_list), region, tail=self.annotation.span_size)

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
            self.logger.error(error_message)
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
        for transcript in self.annotation.get_transcript_assembly(max(span_before, span_after)):
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

    @preconditions(lambda dist: isinstance(dist, int),
                   lambda dist: dist > 0)
    def get_amino_acid_pauses(self, dist, downsample=True):
        """
        Counts the meta-number of 5' mapping positions at the given distance from the specified codon
        Only transcripts with cds of length multiple of 3 are accounted for.
        The only frame in these transcripts is considered.

        :param codon:
        :param span_before:
        :param span_after:
        :return:
        """

        self.logger.info(
            "Counting amino acid specific pauses within 0 to %d nt distance from the first nucleotide of each codon" % dist)
        # amino_acid_count_dict = {}
        amino_acid_count_df = pd.DataFrame(data=0, index=Codons.AMINO_ACID_TABLE.keys(),
                                           columns=range(-1 * dist, 0))

        counter = 1
        wrong_cds_count = 0
        multi_stop_cds_count = 0
        single_stop_cds_count = 0
        no_stop_cds_count = 0

        # FIXME why on earth am I setting a span_size and then suffering to remove them from the ends?
        # FIXME have set all span_sizes to 0. Hope it won't crush
        transcript_assembly = self.annotation.get_transcript_assembly_default_filter(0)
        transcript_count = len(transcript_assembly)
        for transcript in transcript_assembly:
            if counter % np.floor(transcript_count / 10) == 0:
                self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                    counter, floor(100 * (counter - 1) / transcript_count), '%',), )
                self.logger.info("Amount of cds not multiple of 3 is %.2f %s"
                                 % (float(100 * wrong_cds_count) / counter, "%"))
                self.logger.info("Amount of cds with 0, 1, and more STOPs: %d, %d, %d"
                                 % (no_stop_cds_count, single_stop_cds_count, multi_stop_cds_count))
            counter += 1

            count_vector = self.get_count_vector(transcript, span_size=0,
                                                 region=FivePSeqCounts.FULL_LENGTH,
                                                 downsample=downsample)
            # sequence = transcript.get_sequence(self.genome.genome_dict)
            # cds_sequence = sequence[0: len(sequence) - 0]
            # cds_sequence = cds_sequence[transcript.cds_start: transcript.cds_end]
            # count_vector = count_vector[0:len(count_vector) - 0]
            cds_sequence = self.get_cds_sequence_safe(transcript, 0)

            if sum(count_vector) == 0:
                continue

            if len(cds_sequence) != len(count_vector):
                # TODO in the future report the number of problematic transcripts
                # TODO filter these transcripts during plotting
                self.logger.warning("Transcript num %d: cds sequence length %d not equal to count vector length %d"
                                    % (counter, len(cds_sequence), len(count_vector)))
                continue

            if len(count_vector) % 3 != 0:
                wrong_cds_count += 1
                continue

            # NOTE below is the faster implementation based only on non-empty triplets of the vectors
            # NOTE as we don't scroll through the full transcripts, the number of stops will only be counted for non-empty triplets
            # identify 3nt bins with non-zero counts
            ind = np.array(range(0, len(count_vector), 3))
            hits = [sum(count_vector[i:i + 3]) > 0 for i in ind]
            non_empty_ind = ind[hits]

            stop_pos = []
            num_stops = 0
            # loop through non-empty triplets only
            for i in non_empty_ind:
                # loop through all amino acids 20 nucleotides downstream (seven amino-acids)
                for j in range(i + 3, i + 3 + dist, 3):
                    if j + 3 > len(cds_sequence):
                        break
                    codon = cds_sequence[j: j + 3].upper()

                    # NOTE the comparison is case-sensitive and the low-case letters are now not counted
                    # NOTE however, low-case may indicate repetitive regions and it may be advantagous to skip them
                    if (len(codon) == 3) & (codon in Codons.CODON_TABLE.keys()):
                        aa = Codons.CODON_TABLE[codon]
                        if codon in Codons.stop_codons:
                            num_stops += 1
                            if j not in stop_pos:
                                stop_pos.append(j)
                        for p in range(0, 3):
                            d = i - j + p
                            if -1 * d <= dist:
                                amino_acid_count_df.at[aa, d] += count_vector[i + p]

            num_stops = len(stop_pos)
            if num_stops > 1:
                multi_stop_cds_count += 1
            elif num_stops == 1:
                single_stop_cds_count += 1
            else:
                no_stop_cds_count += 1

            # TODO test for problems with the fast implementation and remove the code below if successful
            # comment out the previous implentation
            """
            if len(count_vector) % 3 == 0:
                num_stops = 0
                for i in range(0, len(count_vector), 3):
                    codon = cds_sequence[i: i + 3]

                    # NOTE the comparison is case-sensitive and the low-case letters are now not counted
                    # NOTE however, low-case may indicate repetitive regions and it may be advantagous to skip them
                    if (len(codon) == 3) & (codon in Codons.CODON_TABLE.keys()):
                        aa = Codons.CODON_TABLE[codon]
                        if codon in Codons.stop_codons:
                            num_stops += 1
                        if i > dist:
                            amino_acid_count_df.loc[aa, -1 * dist: 0] += count_vector[i - dist: i]
                        else:
                            amino_acid_count_df.loc[aa, -1 * i: 0] += count_vector[0: i]
                if num_stops > 1:
                    multi_stop_cds_count += 1
                elif num_stops == 1:
                    single_stop_cds_count += 1
                else:
                    no_stop_cds_count += 1
            else:
                wrong_cds_count += 1

        """
        self.logger.debug("Amount of cds not multiple of 3 is %.2f %s"
                          % (float(100 * wrong_cds_count) / counter, "%"))
        self.logger.debug("Amount of cds with 0, 1, and more STOPs: %d, %d, %d"
                          % (no_stop_cds_count, single_stop_cds_count, multi_stop_cds_count))

        # for d in range(1,dist+1):
        #    if (i > d) & (len(cds_sequence) > d):
        #        amino_acid_count_df.loc[aa, -1*d] += count_vector[i - d]
        # amino_acid_count_dict[aa][dist-d-1] += count_vector[i - d]

        return amino_acid_count_df

    @preconditions(lambda dist_from: isinstance(dist_from, int),
                   lambda dist_from: dist_from < 0,
                   lambda dist_to: isinstance(dist_to, int),
                   lambda dist_to: dist_to >= 0)
    def get_codon_pauses(self, dist_from=-30, dist_to=6, downsample=True):
        """
        Counts the meta-number of 5' mapping positions at the given distance from the specified codon
        Only transcripts with cds of length multiple of 3 are accounted for.
        The only frame in these transcripts is considered.

        :param codon:
        :param dist_from: negative distance from each codon
        :param dist_to: positive distance after each codon
        :return:
        """

        self.logger.info(
            "Counting codon specific pauses within %d to %d nt distance from the first nucleotide of each codon" %
            (dist_from, dist_to))

        codon_count_df = pd.DataFrame(data=0, index=Codons.CODON_TABLE.keys(),
                                           columns=range(dist_from, dist_to))

        counter = 1

        transcript_assembly = self.annotation.get_transcript_assembly_default_filter(span_size=dist_to)
        transcript_count = len(transcript_assembly)
        for transcript in transcript_assembly:
            if counter % np.floor(transcript_count / 10) == 0:
                self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                    counter, floor(100 * (counter - 1) / transcript_count), '%',), )
            counter += 1

            count_vector = self.get_count_vector(transcript, span_size=dist_to,
                                                 region=FivePSeqCounts.FULL_LENGTH,
                                                 downsample=downsample)
            count_vector = count_vector[dist_to:len(count_vector)-dist_to]
            cds_sequence = self.get_cds_sequence_safe(transcript, dist_to)

            if sum(count_vector) == 0:
                continue

            if len(cds_sequence) != len(count_vector):
                self.logger.warning("Transcript num %d: cds sequence length %d not equal to count vector length %d"
                                    % (counter, len(cds_sequence), len(count_vector)))
                continue

            # identify 3nt bins with non-zero counts
            ind = np.array(range(0, len(count_vector), 3))
            hits = [sum(count_vector[i:i + 3]) > 0 for i in ind]
            non_empty_ind = ind[hits]

            # loop through non-empty triplets only
            for i in non_empty_ind:
                # loop through all codons dist_from nucleotides downstream and dist_to nucleotides upstream
                j_range = list(np.arange(i, i-dist_to, -3))[::-1] + list(np.arange(i+3, i + 3 - dist_from, 3))
                for j in j_range:
                    if j < 0:
                        continue
                    if j + 3 > len(cds_sequence):
                        break
                    codon = cds_sequence[j: j + 3].upper()

                    if (len(codon) == 3) & (codon in Codons.CODON_TABLE.keys()):
                        for p in range(0, 3):
                            d = i - j + p
                            try:
                                codon_count_df.at[codon, d] += count_vector[i + p]
                            except Exception as e:
                                self.logger.warn("Index out of range: i: %d, j: %d, p: %d, d: %d. %s"
                                                 % (i, j, p, d, str(e)))


        # rename codon_count_df indices by adding amino acid names
        new_index = [Codons.CODON_TABLE.get(codon) + '_' + codon for codon in codon_count_df.index]
        codon_count_df.index = new_index
        codon_count_df = codon_count_df.reindex(sorted(new_index))

        return codon_count_df

    @preconditions(lambda padding: isinstance(padding, int),
                   lambda padding: padding > 0,
                   lambda loci_file: str)
    def get_pauses_from_loci(self, loci_file, padding=100):
        """
        Counts the meta-number of 5' mapping positions at the given distance from the specified loci

        The loci file should contain one locus per row.
        Two tab separated columns should indicate chromosome number and position.

        The distance of 5' mapping positions from each loci is counted within each cds.
        The padding sizes are subtracted from the start and end of each transcript.

        :param padding: int: padding, bp (not to count the first and last regions in the transcripts)
        :param loci_file: str: full path to the file specifying the loci.
        :return:
        """

        self.logger.info(
            "Counting pauses from loci given in file %s" % loci_file)

        loci = pd.read_csv(loci_file, sep="\t", header=None, names=['chr', 'pos'], index_col=None)

        # the results will be kept in a dictionary:
        #   key - distance from any locus
        #   value - number of mapping positions at key distance from any locus
        loci_pauses_dict = {}

        counter = 0
        loci_row = 0

        done = False
        move_transcript = True
        move_locus = False
        tg = self.annotation.get_transcript_assembly_default_filter(0)
        transcript = None
        while (True):
            if counter % 1000 == 0:
                self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                    counter, floor(100 * (counter - 1) / self.annotation.transcript_count), '%',), )

            if move_locus:
                if loci.shape[0] == loci_row:
                    self.logger.debug("Reached the end of loci file (row %d)" % loci_row)
                    break
                loci_row += 1
                move_locus = False
                continue

            if move_transcript:
                try:
                    transcript = tg[counter]
                except:
                    self.logger.debug("Reached the end of transcript assembly (counter: %d)" % counter)
                    break
                counter += 1
                move_transcript = False
                continue

            # check if the locus at the cursor is within the current transcript
            if loci_row < loci.shape[0]:
                if str(transcript.chrom) == str(loci.loc[loci_row, "chr"]):

                    if transcript.cds_genome_start > loci.loc[loci_row, "pos"]:
                        move_locus = True
                        continue

                    elif transcript.cds_genome_end < loci.loc[loci_row, "pos"]:
                        move_transcript = True
                        continue

                    else:
                        count_vector = self.get_count_vector(transcript, 0,
                                                             FivePSeqCounts.FULL_LENGTH,
                                                             downsample=True)

                        if len(count_vector) > 2 * padding:
                            for i in range(padding, len(count_vector) - padding):
                                if count_vector[i] > 0:
                                    genomic_position = transcript.cds_genome_start + i
                                    distance = genomic_position - loci.loc[loci_row, "pos"]
                                    if distance in loci_pauses_dict.keys():
                                        loci_pauses_dict[distance] += count_vector[i]
                                    else:
                                        loci_pauses_dict.update({distance: count_vector[i]})
                            move_locus = True
                        else:
                            move_locus = True

                elif str(transcript.chrom) > str(loci.loc[loci_row, "chr"]):
                    move_locus = True
                    continue

                else:
                    move_transcript = True
                    continue
            else:
                break
        # turn the dictionary into a metacount vector, with indices from -1*maxdistance to 0
        self.logger.debug("Merging the dictionary into metacounts")
        maxdist = max(loci_pauses_dict.keys())
        metacount_vector = [0] * maxdist
        for i in range(-1 * maxdist, 0):
            if i in loci_pauses_dict.keys():
                metacount_vector[i] = loci_pauses_dict[i]
        metacount_series = pd.Series(data=metacount_vector, index=np.arange(-1 * maxdist, 0))

        return metacount_series

    @preconditions(lambda num: isinstance(num, int))
    def top_populated_transcript_indices(self, num=1000):
        """
        Returns indices of top populated transcripts.
        A populated transcript is defined as the one with most length-relative number
        of positions with non-zero counts.

        :param num: int: number of transcript indices to return
        :return: [int]: a list of transcript indices in the transcript assembly
        """
        populated = [0] * self.annotation.transcript_count
        for i in range(self.annotation.transcript_count):
            transcript = self.annotation.transcript_assembly[i]
            count_vector = self.get_count_vector(transcript, 0, FivePSeqCounts.FULL_LENGTH, downsample=False)
            populated[i] = sum(count_vector > 0) / len(count_vector)
        populated_indices = sorted(range(len(populated)), key=lambda k: populated[k])

        return populated_indices


class CountManager:
    """
    This module implements a set of static functions to handle count vectors retrieved from FivePSeqCounts class.
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
        Takes a vector of position-wise int counts across full length transcripts
        and returns counts for three different frames from 0 to 2,
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

        elif region == FivePSeqCounts.TERM:
            # NOTE the frames relative to START and TERM should be aligned in the future
            # NOTE (if cds length is not a multiple of 3)
            frame0_array = [count_vector[i] for i in list(reversed(range(len(count_vector) - 3 - tail, -1 + tail, -3)))]
            frame1_array = [count_vector[i] for i in list(reversed(range(len(count_vector) - 2 - tail, -1 + tail, -3)))]
            frame2_array = [count_vector[i] for i in list(reversed(range(len(count_vector) - 1 - tail, -1 + tail, -3)))]

        else:
            error_message = "Invalid region %s specified: should be either %s or %s" \
                            % (region, FivePSeqCounts.START, FivePSeqCounts.TERM)
            logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)
            logger.error(error_message)
            raise Exception(error_message)

        return frame0_array, frame1_array, frame2_array

    @staticmethod
    @preconditions(lambda count_vector: isinstance(count_vector, list),
                   lambda count_vector: isinstance(count_vector[0], int),
                   lambda region: isinstance(region, str),
                   lambda tail: isinstance(tail, int),
                   lambda tail: tail >= 0)
    def count_vector_to_series(count_vector, region, tail=0):
        """
        Takes a vector of counts, indexes them with distances from the specified region.
        Returns a series with indexes as genomic coordinates from start/stop codons, and values as counts at each coordinates.

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

        :param count_vector: [int]: a vector of summed position-wise counts for a (meta)-transcript
        :param region: str: the region respective to which the distance is calculated
        :param tail: int:

        :return: pandas.Series: a series with indices as genomic coordinates* and values as meta counts.
        *-corresponding to positions' distance from nucleotide 0 of START/TERM codons

        """

        if region == FivePSeqCounts.START:
            d = np.arange(-tail, len(count_vector) - tail)

        elif region == FivePSeqCounts.TERM:
            d = np.arange(-(len(count_vector) - tail - 3), tail + 3)

        else:
            error_message = "Invalid region %s specified: should be either %s or %s" \
                            % (region, FivePSeqCounts.START, FivePSeqCounts.TERM)
            logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)
            logger.error(error_message)
            raise Exception(error_message)

        counts_series = pd.Series(data=count_vector, index=d)

        return counts_series

    @staticmethod
    @preconditions(lambda count_vector: isinstance(count_vector, list),
                   lambda count_vector: isinstance(count_vector[0], int),
                   lambda region: isinstance(region, str),
                   lambda tail: isinstance(tail, int),
                   lambda tail: tail >= 0)
    def count_vector_to_df(count_vector, region, tail=0):
        """
        Takes a vector of counts, indexes them with distances from the specified region.
        Returns a dataframe with indexes as genomic coordinates from start/stop codons, and values as counts at each coordinates.

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

        :param count_vector: [int]: a vector of summed position-wise counts for a (meta)-transcript
        :param region: str: the region respective to which the distance is calculated
        :param tail: int:

        :return: pandas.Series: a series with indices as genomic coordinates* and values as meta counts.
        *-corresponding to positions' distance from nucleotide 0 of START/TERM codons

        """

        if region == FivePSeqCounts.START:
            d = np.arange(-tail, len(count_vector) - tail)

        elif region == FivePSeqCounts.TERM:
            d = np.arange(-(len(count_vector) - tail - 3), tail + 3)

        else:
            error_message = "Invalid region %s specified: should be either %s or %s" \
                            % (region, FivePSeqCounts.START, FivePSeqCounts.TERM)
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error(error_message)
            raise Exception(error_message)

        counts_df = pd.DataFrame({'D': d, 'C': count_vector})

        return counts_df

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

        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug(
            "Retrieving count-sums per frame relative to %s ..." % region)

        # Create an empty dataframe
        n = len(count_vector_list)
        frame_counts_df = pd.DataFrame({'F0': [0] * n,
                                        'F1': [0] * n,
                                        'F2': [0] * n})
        for t_ind in range(0, n):
            # Print status update to console
            if t_ind % 100 == 0:
                print "\r>>Transcript count: %d (%d%s)\t" % (
                    t_ind, np.floor(100 * (t_ind - 1) / n), '%'),

            # extract frame count vectors from count vectors
            count_vector = count_vector_list[t_ind]
            frame_counts = CountManager.extract_frame_count_vectors(count_vector, span_size, region)

            # sum-up counts in each frame and add to the dataframe
            for f in range(0, 3):
                frame_counts_df.iloc[t_ind, f] = sum(frame_counts[f])

        return frame_counts_df

    @staticmethod
    @preconditions(lambda file_path: isinstance(file_path, str))
    def read_counts_as_list(file_path):
        """
        Reads and returns a list of count vectors, each corresponding to a transcript.

        :param file_path: str: full path to the file
        :return: [[int]]: list of count vectors (a count vector is a list of int counts)
        :exception: raises IOError if file does not exist
        """

        if not os.path.exists(file_path):
            error_message = "Problem reading counts: the file %s does not exist" % file_path
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug("Reading count file %s" % file_path)
        df = pd.read_csv(file_path, header=None, sep="|")
        count_vector_list = [[]] * len(df)
        for i in range(0, len(df)):
            count_vector_list[i] = map(int, df.iloc[i, 0].split("\t"))
        return count_vector_list

    @staticmethod
    @preconditions(lambda file: isinstance(file, str))
    def read_meta_counts(file):
        """
        Reads the meta count file as pandas DataFrame.
        These files are saved with tab separator.
        They have two columns, but no column names.
        This function assigns names to read DataFrame:
            D (distance from START/TERM)
            C (meta counts)
        The number of rows corresponds to 2*span_size

        :param file: str: full path to the file
        :return: pandas DataFrame: a dataframe with D and C columns
        :exception: raises IOError if file does not exist
        """
        if not os.path.exists(file):
            error_message = "Problem reading meta counts: the file %s does not exist" % file
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug("Reading meta count file %s" % file)
        meta_count = pd.read_csv(file, sep="\t", header=None, names=["D", "C"])
        return meta_count

    @staticmethod
    @preconditions(lambda file: isinstance(file, str))
    def read_frame_counts(file):
        """
        Reads the frame count file as pandas DataFrame.
        The file has a header with four columns:
        (no name: transcript number), F0, F1, F2
        A four-column dataFrame is created and returned accordingly.
        The number of rows corresponds to the number of transcripts.

        :param file: str: full path to the file
        :return: pandas DataFrame: a dataframe with transcript number and F0, F1, F2 frame counts
        :exception: raises IOError if file doesn't exist
        """
        if not os.path.exists(file):
            error_message = "Problem reading frame counts: the file %s does not exist" % file
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug("Reading frame counts file %s" % file)

        frame_counts = pd.read_csv(file, sep="\t")
        return frame_counts

    @staticmethod
    @preconditions(lambda file: isinstance(file, str))
    def read_amino_acid_df(file):
        """
        Reads a pandas DataFrame from amino acid pauses file.
        The file is stored with a header indicating distance from amino acids.
        The file has row names indicating names of amino acids.
        The dataFrame is read with indicated columns and row names.

        :param file: str: full path to file
        :return: pandas DataFrame: index is amino acids, columns - distance
        :exception: raises IOError if file does not exist
        """
        if not os.path.exists(file):
            error_message = "Problem reading amino acid pauses: the file %s does not exist" % file
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug("Reading amino acids pauses file %s" % file)

        amino_acid_df = pd.read_csv(file, sep="\t", header=0, index_col=0)

        return amino_acid_df

    @staticmethod
    def top_populated_count_vector_indices(count_vector_list, num=1000):
        """
        Returns indices of top populated count_vectors (transcripts).
        A populated count_vector (transcript) is defined as the one with most length-relative number
        of positions with non-zero counts.

        :param count_vector_list: [[int]]: a list of count vectors
        :param num: int: number of indices to return
        :return: [int]: a list of count_vector indices in the count_vector_list
        """

        populated = [0] * len(count_vector_list)
        for i in range(len(count_vector_list)):
            count_vector = count_vector_list[i]
            populated[i] = float(np.count_nonzero(count_vector)) / len(count_vector)
        populated_indices = sorted(range(len(populated)), key=lambda k: populated[k], reverse=True)

        return populated_indices[0:num]

    @staticmethod
    def canonical_transcript_indices(count_dir):
        """
        Reads and returns canonical transcript indices from the canonical transcript indices file,
        if such a file exists.

        :return: [int] indices of transcript with canonical start and stop codons or None if no such file exists.

        """

        canonical_index_file = os.path.join(count_dir, FivePSeqOut.CANONICAL_TRANSCRIPT_INDEX_FILE)
        if os.path.exists(canonical_index_file):
            transcript_index = list(pd.read_csv(canonical_index_file, header=None).iloc[:, 0])
            return transcript_index
        else:
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug(
                "Problem retrieving canonical transcript indices. No file %s exists. "
                "The filter will return None." % canonical_index_file)
            return None

    @staticmethod
    @preconditions(lambda file_path: isinstance(file_path, str))
    def read_count_vector(file_path):
        """
        Reads and returns a list of counts from a new-line separated file of counts.

        :param file_path: str: full path to the file
        :return: [int]: a list of counts
        :exception: raises IOError if file does not exist
        """

        if not os.path.exists(file_path):
            error_message = "Problem reading counts: the file %s does not exist" % file_path
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).debug("Reading count file %s" % file_path)
        if os.stat(file_path).st_size == 0:
            counts = []
        else:
            counts = pd.read_csv(file_path, header=None)
            counts = map(int, counts.iloc[:, 0])

        return counts
