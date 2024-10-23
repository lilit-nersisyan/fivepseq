import collections
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

    START = "START"
    TERM = "STOP"
    MID = "MIDDLE"
    FULL_LENGTH = "full_length"
    ALL = "all"

    START_CODON = "start"
    STOP_CODON = "stop"
    TRANSCRIPT_LENGTH = "len"
    TRANSCRIPT_3NT = "3nt"
    NUMBER_READS = "NumOfReads"
    NUMBER_READS_DOWNSAMPLED = "NumOfReadsDownsampled"
    NUMBER_POSITIONS = "NumOfMapPositions"

    COUNT_THRESHOLD = 100
    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

    count_distribution_dict = None
    outlier_lower = None
    downsample_constant = None
    outlier_probability = None

    config = None
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
    codon_genome_usage_df = None
    codon_count_df = None
    amino_acid_count_df = None
    amino_acid_island_count_df = None
    dicodon_count_df = None
    dipeptide_count_df = None
    tricodon_count_df = None
    tripeptide_count_df = None
    codon_stats_df = None
    amino_acid_stats_df = None

    amino_acid_genome_usage_df = None

    start_codon_dict = None
    stop_codon_dict = None
    canonical_transcript_index = None
    transcript_descriptors = None
    outliers = None

    is_geneset = False

    loci_overlaps = None
    READ_LOCATIONS_ALL = "_ALL"
    READ_LOCATIONS_3UTR = "_3UTR"
    READ_LOCATIONS_5UTR = "_5UTR"
    READ_LOCATIONS_CDS = "_CDS"

    MASK_DIST = 20
    TRIPEPTIDE_POS = -11
    DIPEPTIDE_POS = -14
    CSPAN = -30

    missing_chroms = []

    def __init__(self, alignment, annotation, genome, config, downsample_constant, is_geneset=False,
                 transcript_filter=None):
        """
        Initializes a FivePSeqCounts object with Alignment and Annotation instances.

        :param alignment: fivepseq.logic.structures.Alignment type object
        :param annotation: fivepseq.logic.structures.Annotation type object
        :param genome: fivepseq.logic.structures.Genome: Genome type object
        :param outlier_probability: a float setting the probability threshold for Poisson distribution that will be used to downsample outliers
        :param downsample_constant: a float specifying a constant threshold: higher values will be down-sampled to this constant (without Poisson check)
        """

        self.alignment = alignment
        self.annotation = annotation
        self.genome = genome
        self.transcript_filter = transcript_filter
        self.config = config
        self.outlier_probability = config.args.op
        self.outlier_lower = downsample_constant
        self.outliers = []
        self.start_codon_dict = {}
        self.stop_codon_dict = {}
        self.canonical_transcript_index = []
        self.is_geneset = is_geneset
        self.fpat = None

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

        transcript_assembly = self.annotation.get_transcript_assembly(span_size=0)
        transcript_count = len(transcript_assembly)
        self.transcript_descriptors = pd.DataFrame(data=None,
                                                   index=range(transcript_count),
                                                   columns=[self.START_CODON,
                                                            self.STOP_CODON,
                                                            self.TRANSCRIPT_LENGTH,
                                                            self.TRANSCRIPT_3NT,
                                                            self.NUMBER_READS,
                                                            self.NUMBER_READS_DOWNSAMPLED,
                                                            self.NUMBER_POSITIONS])

        count_distribution_dict = {}

        for transcript_ind in range(transcript_count):
            transcript = transcript_assembly[transcript_ind]

            cds_sequence = self.get_cds_sequence_safe(transcript, 0)
            count_vector = self.get_count_vector_safe(transcript, 0)

            # NOTE the count distribution does not include values 0 to avoid skewness for outlier detection
            for c in count_vector:
                if c > 0:
                    if c in count_distribution_dict:
                        count_distribution_dict[c] += 1
                    else:
                        count_distribution_dict[c] = 1

            start_codon = cds_sequence[0:3]
            stop_codon = cds_sequence[len(cds_sequence) - 3:len(cds_sequence)]

            if (start_codon == codons.Codons.START_CODON) & (stop_codon in codons.Codons.stop_codons):
                self.canonical_transcript_index.append(transcript_ind)

            self.transcript_descriptors.at[transcript_ind, self.START_CODON] = start_codon
            self.transcript_descriptors.at[transcript_ind, self.STOP_CODON] = stop_codon
            self.transcript_descriptors.at[transcript_ind, self.TRANSCRIPT_3NT] = str(len(cds_sequence) % 3 == 0)
            self.transcript_descriptors.at[transcript_ind, self.TRANSCRIPT_LENGTH] = len(cds_sequence)
            self.transcript_descriptors.at[transcript_ind, self.NUMBER_READS] = int(np.sum(count_vector))
            self.transcript_descriptors.at[transcript_ind, self.NUMBER_POSITIONS] = np.count_nonzero(count_vector)

            if start_codon in self.start_codon_dict.keys():
                self.start_codon_dict[start_codon] += 1
            else:
                self.start_codon_dict.update({start_codon: 1})

            if stop_codon in self.stop_codon_dict.keys():
                self.stop_codon_dict[stop_codon] += 1
            else:
                self.stop_codon_dict.update({stop_codon: 1})

        self.count_distribution_dict = collections.OrderedDict(sorted(count_distribution_dict.items()))
        self.outlier_lower = self.get_outlier_lower()

        self.logger.info("The lower bound for outliers set as %f " % self.outlier_lower)

        # also store downsampled transcript counts
        for transcript_ind in range(transcript_count):
            transcript = transcript_assembly[transcript_ind]

            count_vector_downsampled = self.get_count_vector(transcript, span_size=0,
                                                             region=self.FULL_LENGTH, downsample=True)
            self.transcript_descriptors.at[transcript_ind, self.NUMBER_READS_DOWNSAMPLED] = int(
                np.sum(count_vector_downsampled))

        self.logger.info("Done generating transcript descriptors")

    def get_count_distribution_dict(self):
        return self.count_distribution_dict

    def get_count_distribution(self):

        if self.count_distribution_dict is None:
            self.generate_transcript_descriptors()

        count_distribution = []
        for c, f in self.count_distribution_dict.items():
            for i in range(f):
                count_distribution.append(c)

        return count_distribution

    def set_count_distribution_dict(self, count_distribution_dict):
        """
        Sets the count distribution according to the specified count vector.

        :param count_distribution_dict: an ordered dictionary of count frequencies
        :return:
        """
        if len(count_distribution_dict) == 0:
            self.count_distribution_dict = None
        else:
            self.count_distribution_dict = count_distribution_dict

    def get_outlier_lower(self):
        """
        Returns the lower bound for outliers detected as points lying self.downsample_by number times higher than the
        25-75% interquartile range.

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
        ind = np.where(np.asarray(ps) <= self.outlier_probability)[0].tolist()

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

        logging.getLogger(config.FIVEPSEQ_LOGGER).info("Generating count vectors")

        if self.count_vector_list_full_length is not None:
            if self.count_vector_list_term is not None:
                if self.count_vector_list_start is not None:
                    logging.getLogger(config.FIVEPSEQ_LOGGER).warning("All count vectors are already generated")

        # otherwise, retrieve the counts from the alignment file, referencing the transcript assembly
        self.logger.info("Retrieving counts (span size :%d)..."
                         % self.annotation.span_size)

        # initialize empty vectors
        transcript_count = len(self.annotation.get_transcript_assembly())
        self.count_vector_list_full_length = [None] * transcript_count
        self.count_vector_list_term = [None] * transcript_count
        self.count_vector_list_start = [None] * transcript_count

        # setup the the counter
        counter = 1
        ta = self.annotation.get_transcript_assembly()
        for i in range(transcript_count):
            transcript = ta[i]

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

    @preconditions(lambda span_size: isinstance(span_size, int),
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
            if downsample and any(x > self.outlier_lower for x in count_vector):
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
        count_vector = list(map(int, count_vector))

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

                logging.getLogger(config.FIVEPSEQ_LOGGER). \
                    debug("Transcript %s at the beginning of the genome padded with %d zeros"
                          % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))

            else:
                t_len = transcript.spanning_segment.end - transcript.spanning_segment.start
                diff = transcript.spanning_segment.end - len(self.genome.genome_dict[transcript.chrom].seq)

                if diff > span_size:
                    # NOTE wrongly annotated transcripts go outside genome boundaries,
                    # NOTE return an empty vector spanned by span size as a safe way of discarding such transcripts
                    count_vector = [0] * t_len
                    logging.getLogger(config.FIVEPSEQ_LOGGER). \
                        debug("Transcript %s exceeds genome dimensions by %d bases"
                              % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))

                else:
                    t_subchain = transcript.get_subchain(diff, transcript.spanning_segment.end, stranded=False)

                    subchain_counts = list(t_subchain.get_counts(self.alignment.bam_array))
                    count_vector = subchain_counts + [0] * diff

                    logging.getLogger(config.FIVEPSEQ_LOGGER). \
                        debug("Transcript %s at the end of the genome padded with %d zeros"
                              % (FivePSeqOut.get_transcript_attr(transcript, "Name"), diff))


        count_vector = self.filter_pattern(transcript, count_vector, span_size)

        return count_vector

    def get_sequence(self, transcript, transcript_span_size, desired_span_size):
        if desired_span_size > transcript_span_size:
            raise ValueError("Desired span size %d bigger than the transcript span size %d"
                             % (desired_span_size, transcript_span_size))

        try:
            sequence = transcript.get_sequence(self.genome.genome_dict)
            desired_seq = sequence[transcript.cds_start + transcript_span_size - desired_span_size:
                                   transcript.cds_end + transcript_span_size + desired_span_size]
        except:
            t_len = transcript.cds_end - transcript.cds_start
            desired_seq = ''.join(['N'] * t_len + 2 * desired_span_size)

        return desired_seq

    def get_cds_sequence_safe(self, transcript, span_size):
        # NOTE a dangerous code here. Works correctly only if the input span size is the same as in the transcript.
        # TOCHANGE
        try:
            sequence = transcript.get_sequence(self.genome.genome_dict)
            cds_sequence = sequence[transcript.cds_start + span_size: transcript.cds_end + span_size]
        except:
            if transcript.chrom not in self.genome.genome_dict.keys():
                if transcript.chrom not in self.missing_chroms:
                    self.missing_chroms.append(transcript.chrom)
                    logging.getLogger(config.FIVEPSEQ_LOGGER).warn(
                        "No chromosome named %s found in the genome sequence" % transcript.chrom)
                t_len = transcript.spanning_segment.end - transcript.spanning_segment.start
                cds_sequence = ''.join(['N'] * t_len)

            elif transcript.spanning_segment.start < 0:
                diff = -1 * transcript.spanning_segment.start
                t_subchain = transcript.get_subchain(diff, transcript.spanning_segment.end, stranded=False)
                sequence = t_subchain.get_sequence(self.genome.genome_dict)

                if span_size < diff:
                    cds_sequence = sequence[transcript.cds_start + span_size - diff: transcript.cds_end + span_size]
                else:  # TODO I don't know how to get sequence in this case: need debugging
                    cds_sequence = sequence[transcript.cds_start + span_size - diff: transcript.cds_end + span_size]

                logging.getLogger(config.FIVEPSEQ_LOGGER). \
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

    def get_amino_acid_pauses(self):
        if self.amino_acid_count_df is None:
            self.compute_codon_pauses(dist_from=self.get_cspan())

        return self.amino_acid_count_df

    def get_amino_acid_island_pauses(self):
        if self.amino_acid_island_count_df is None:
            self.compute_aa_islands(dist_from=self.get_cspan())

        return self.amino_acid_island_count_df

    def get_codon_pauses(self):
        if self.codon_count_df is None:
            self.compute_codon_pauses(dist_from=self.get_cspan())

        return self.codon_count_df

    def get_tricodon_pauses(self):
        if self.tricodon_count_df is None:
            self.compute_codon_pauses(dist_from=self.get_cspan())

        return self.tricodon_count_df

    def get_dicodon_pauses(self):
        if self.dicodon_count_df is None:
            self.compute_codon_pauses(dist_from=self.get_cspan())

        return self.dicodon_count_df

    def get_dipeptide_pauses(self):
        if self.dipeptide_count_df is None:
            self.compute_codon_pauses(dist_from=self.get_cspan())

        return self.dipeptide_count_df

    def get_tripeptide_pauses(self):
        if self.tripeptide_count_df is None:
            self.compute_codon_pauses(dist_from=self.get_cspan())

        return self.tripeptide_count_df

    def store_out_of_frame_stats(self, period, asite_frag=None, downsample=True):
        """
        Stores the cummulative counts located A-site distance from each codon in-frame and out-of-frame.

        :param: period: int a number indicating the length of the periodicity to search for
        :param: asite_frag: int length from A site to ribosome protection site (computed from period if not provided)
        :return: df: pd.DataFrame: a dataframe containing columns ["f0", "f1", "f2"] and rows for each codon
        """

        logging.info("Counting codon pauses in- and out-of- frame ")

        if asite_frag is None:
            asite_frag = period - 13  # inferred from 30 and -17 for yeast, but can be changed for other species

        logging.info("Options: A-site distance = %d" % asite_frag)

        out_of_frame_df = pd.DataFrame(data=0, index= Codons.CODON_TABLE.keys(),
                                       columns=["f0", "f1", "f2"])


        span_size = self.annotation.span_size
        cnt = 0
        ta = self.annotation.get_transcript_assembly(span_size=span_size)
        transcript_count = len(ta)
        for t in range(transcript_count):
            transcript = ta[t]

            if np.floor(transcript_count / 1000) > 0 and t % 1000 == 0:
                logging.info("\r>>Transcript count: %d (%d%s)\t" % (
                    t, floor(100 * (t - 1) / transcript_count), '%',), )

            count_vector = self.get_count_vector(transcript, span_size=0,
                                                 region=FivePSeqCounts.FULL_LENGTH,
                                                 downsample=downsample)

            cds_sequence = self.get_cds_sequence_safe(transcript, 0)

            if sum(count_vector) == 0:
                continue

            if len(cds_sequence) != len(count_vector):
                self.logger.warning("Transcript num %d: cds sequence length %d not equal to count vector length %d"
                                    % (t, len(cds_sequence), len(count_vector)))
                continue

            for f in [0, 1, 2]:

                # loop through codons in the frame
                for i in range(f, len(count_vector), 3):
                    codonA = cds_sequence[i: i + 3].upper()

                    if not codonA in out_of_frame_df.index:
                        break

                    if i < asite_frag:
                        continue

                    out_of_frame_df.loc[codonA, "f" + str(f)] += count_vector[i - asite_frag]

        # rename out_of_frame_df indices by adding amino acid names
        new_index = [Codons.CODON_TABLE.get(codon) + '_' + codon for codon in out_of_frame_df.index]
        out_of_frame_df.index = new_index
        out_of_frame_df = out_of_frame_df


        return out_of_frame_df

    # TODO a modification for di-codon counts to me incorporated inseat of get_codon_pauses in the future
    @preconditions(lambda dist_from: isinstance(dist_from, int),
                   lambda dist_from: dist_from < 0,
                   lambda dist_from: -1*dist_from % 3 == 0,
                   lambda dist_to: isinstance(dist_to, int),
                   lambda dist_to: dist_to >= 0)
    def compute_codon_pauses(self, dist_from=-30, dist_to=9, downsample=True):
        """
        Counts the meta-number of 5' mapping positions at the given distance from a codon or codon-pair
        Only transcripts with cds of length multiple of 3 are accounted for.
        The only frame in these transcripts is considered.

        :param codon:
        :param dist_from: negative distance from each codon or codon-pair
        :param dist_to: positive distance after each codon or codon-pair
        :param mask_dist: the number of positions to mask in the beginning and end of the gene body

        :return:
        """

        self.logger.info(
            "Counting codon specific pauses within %d to %d nt distance from the first nucleotide of each codon" %
            (dist_from, dist_to))

        mask_dist = self.get_mask_dist()

        codon_count_df = pd.DataFrame(data=0, index=Codons.CODON_TABLE.keys(),
                                      columns=range(dist_from, dist_to))

        dicodon_count_df = pd.DataFrame(data=0, index=Codons.get_dicodon_table().keys(),
                                        columns=range(dist_from + 3, dist_to + 3))

        dipeptide_count_df = pd.DataFrame(data=0, index=Codons.get_dipeptide_list(),
                                          columns=range(dist_from + 3, dist_to + 3))

        tricodon_count_df = pd.DataFrame(data=0, index=Codons.get_tricodon_table().keys(),
                                         columns=range(dist_from + 6, dist_to + 6))

        tripeptide_count_df = pd.DataFrame(data=0, index=Codons.get_tripeptide_list(),
                                           columns=range(dist_from + 6, dist_to + 6))

        self.codon_genome_usage_df = pd.DataFrame(data=0, index=Codons.CODON_TABLE.keys(),
                                                  columns=['abs', 'fraction'])
        self.amino_acid_genome_usage_df = pd.DataFrame(data=0, index=Codons.AMINO_ACID_TABLE.keys(),
                                                       columns=['abs', 'fraction'])

        counter = 1

        transcript_assembly = self.annotation.get_transcript_assembly(
            span_size=0)  # don't take more than the gene body (different from previous versions)
        transcript_count = len(transcript_assembly)
        for t in range(transcript_count):
            transcript = transcript_assembly[t]
            if np.floor(transcript_count / 1000) > 0 and counter % 1000 == 0:
                self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                    counter, floor(100 * (counter - 1) / transcript_count), '%',), )
            counter += 1

            count_vector = self.get_count_vector(transcript, span_size=0,
                                                 region=FivePSeqCounts.FULL_LENGTH,
                                                 downsample=downsample)

            cds_sequence = self.get_cds_sequence_safe(transcript, 0)

            if sum(count_vector) == 0:
                continue

            if len(cds_sequence) != len(count_vector):
                self.logger.warning("Transcript num %d: cds sequence length %d not equal to count vector length %d"
                                    % (counter, len(cds_sequence), len(count_vector)))
                continue

            if (mask_dist >= 3):
                # v1.0b3 mask the first and last 20 counts to avoid initiation affecting codon-specific counts
                count_vector[0:mask_dist] = [0] * mask_dist
                # v1.0b3 mask the last 20 nucleotides to avoid termination affecting codon-specific counts, but keep STOP codon counts
                cds_sequence = cds_sequence[0:len(cds_sequence) - mask_dist] + ''.join(
                    'N' * (mask_dist - 3)) + cds_sequence[len(cds_sequence) - 3:len(cds_sequence)]

            # v1.0b3 add stretches of 0's to count_vector and N's to cds_sequence to avoid checking vector boundaries
            count_vector = [0] * (-1 * dist_from) + count_vector + [0] * dist_to
            cds_sequence = ''.join('N' * (-1 * dist_from)) + cds_sequence + ''.join('N' * dist_to)

            # store genome usage stats
            for i in range(0, len(cds_sequence), 3):
                codon = cds_sequence[i: i + 3].upper()
                if codon in self.codon_genome_usage_df.index:
                    self.codon_genome_usage_df.at[codon, "abs"] += 1
                    amino_acid = Codons.CODON_TABLE.get(codon)
                    self.amino_acid_genome_usage_df.at[amino_acid, "abs"] += 1

            # identify 3nt bins with non-zero counts
            ind = np.array(range(0, len(count_vector), 3))
            hits = [sum(count_vector[i:i + 3]) > 0 for i in ind]
            non_empty_ind = ind[hits]

            # loop through non-empty triplets only
            for i in non_empty_ind:
                # loop through all codons dist_from nucleotides downstream and dist_to nucleotides upstream
                j_range = list(np.arange(i, i - dist_to, -3))[::-1] + list(np.arange(i + 3, i + 3 - dist_from, 3))
                for j in j_range:
                    if j < 0:
                        continue
                    if j + 3 > len(cds_sequence):
                        break
                    codonA = cds_sequence[j: j + 3].upper()

                    if j - 3 >= 0:
                        codonP = cds_sequence[j - 3: j].upper()
                    else:
                        codonP = 'NNN'
                    if j - 6 >= 0:
                        codonE = cds_sequence[j - 6: j - 3].upper()
                    else:
                        codonE = 'NNN'

                    if (len(codonA) == 3) & (codonA in Codons.CODON_TABLE.keys()):
                        for p in range(0, 3):
                            d = i - j + p
                            try:
                                codon_count_df.at[codonA, d] += count_vector[i + p]
                                if len(codonP) == 3 and codonP in Codons.CODON_TABLE.keys():
                                    dicodon_count_df.at[codonP + codonA, d + 3] += count_vector[i + p]
                                    dipeptide = Codons.get_peptide_from_codon_list([codonP, codonA])
                                    dipeptide_count_df.at[dipeptide, d + 3] += count_vector[i + p]

                                    if len(codonE) == 3 and codonE in Codons.CODON_TABLE:
                                        tricodon_count_df.at[codonE + codonP + codonA, d + 6] += count_vector[i + p]
                                        tripeptide = Codons.get_peptide_from_codon_list([codonE, codonP, codonA])
                                        tripeptide_count_df.at[tripeptide, d + 6] += count_vector[i + p]

                            except Exception as e:
                                self.logger.warn("Index out of range: i: %d, j: %d, p: %d, d: %d. %s"
                                                 % (i, j, p, d, str(e)))

        self.codon_genome_usage_df.loc[:, "fraction"] = self.codon_genome_usage_df.loc[:, "abs"] / sum(
            self.codon_genome_usage_df.loc[:, "abs"])
        self.amino_acid_genome_usage_df.loc[:, "fraction"] = self.amino_acid_genome_usage_df.loc[:, "abs"] / sum(
            self.amino_acid_genome_usage_df.loc[:, "abs"])
        self.amino_acid_count_df = self.codon_to_amino_acid_count_df(codon_count_df)
        self.tripeptide_count_df = self.filter_codon_counts(tripeptide_count_df, self.get_tripeptide_pos())
        self.dipeptide_count_df = self.filter_codon_counts(dipeptide_count_df, self.get_dipeptide_pos())

        # rename codon_count_df indices by adding amino acid names
        new_index = [Codons.CODON_TABLE.get(codon) + '_' + codon for codon in codon_count_df.index]
        codon_count_df.index = new_index
        self.codon_count_df = codon_count_df

        # rename codon_count_df indices by adding amino acid names
        self.logger.info("Mapping tricodons to amino acid names")
        tricodon_count_df.index = Codons.get_tricodon_full_index()
        self.tricodon_count_df = self.filter_codon_counts(tricodon_count_df, self.get_tripeptide_pos())

        # rename codon_count_df indices by adding amino acid names
        self.logger.info("Mapping dicodons to amino acid names")
        dicodon_count_df.index = Codons.get_dicodon_full_index()
        self.dicodon_count_df = self.filter_codon_counts(dicodon_count_df, self.get_dipeptide_pos())

        return

    def get_mask_dist(self):
        if self.config.args.no_mask:
            mask_dist = 0
            self.logger.info("Transcript boundaries will not be masked")
        else:
            if hasattr(config.args, "codon_mask_size"):
                mask_dist = config.args.codon_mask_size
            else:
                mask_dist = self.MASK_DIST
            self.logger.info("Transcript boundaries will be masked by %d nucleotides" % mask_dist)
        return mask_dist

    @preconditions(lambda dist_from: isinstance(dist_from, int),
                   lambda dist_from: dist_from < 0,
                   lambda dist_to: isinstance(dist_to, int),
                   lambda dist_to: dist_to >= 0)
    def compute_aa_islands(self, dist_from=-30, dist_to=3, downsample=True):
        """
        Counts the meta-number of 5' mapping positions at the given distance from each amino acid encoding codon
        Only transcripts with cds of length multiple of 3 are accounted for.
        The only frame in these transcripts is considered.

        :param codon:
        :param dist_from: negative distance from each codon
        :param dist_to: positive distance after each codon
        :param mask_dist: the number of positions to mask in the beginning and end of the gene body

        :return:
        """
        # expand dist_to
        dist_to = -1*dist_from

        self.logger.info(
            "Counting amino acid specific pauses at amino acid islands within %d to %d nt distance from the first nucleotide of each codon" %
            (dist_from, dist_to))

        mask_dist = self.get_mask_dist()

        aa_island_count_df = pd.DataFrame(data=0, index=Codons.AMINO_ACID_TABLE.keys(),
                                          columns=range(dist_from, dist_to))

        for aa in Codons.AMINO_ACID_TABLE.keys():
            self.logger.info( "%s islands" % aa )
            aa_codons = Codons.AMINO_ACID_TABLE[aa]

            counter = 1

            transcript_assembly = self.annotation.get_transcript_assembly(
                span_size=0)  # don't take more than the gene body (different from previous versions)
            transcript_count = len(transcript_assembly)
            for t in range(transcript_count):
                transcript = transcript_assembly[t]
                if np.floor(transcript_count / 1000) > 0 and counter % 1000 == 0:
                    self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                        counter, floor(100 * (counter - 1) / transcript_count), '%',), )
                counter += 1

                count_vector = self.get_count_vector(transcript, span_size=0,
                                                     region=FivePSeqCounts.FULL_LENGTH,
                                                     downsample=downsample)

                cds_sequence = self.get_cds_sequence_safe(transcript, 0)

                if sum(count_vector) == 0:
                    continue

                if len(cds_sequence) != len(count_vector):
                    self.logger.warning("Transcript num %d: cds sequence length %d not equal to count vector length %d"
                                        % (counter, len(cds_sequence), len(count_vector)))
                    continue

                if (mask_dist >= 3):
                    # v1.0b3 mask the first and last 20 counts to avoid initiation affecting codon-specific counts
                    count_vector[0:mask_dist] = [0] * mask_dist
                    # v1.0b3 mask the last 20 nucleotides to avoid termination affecting codon-specific counts, but keep STOP codon counts
                    cds_sequence = cds_sequence[0:len(cds_sequence) - mask_dist] + ''.join(
                        'N' * (mask_dist - 3)) + cds_sequence[len(cds_sequence) - 3:len(cds_sequence)]


                # v1.0b3 add stretches of 0's to count_vector and N's to cds_sequence to avoid checking vector boundaries
                count_vector = [0] * (-1 * dist_from) + count_vector + [0] * dist_to
                cds_sequence = ''.join('N' * (-1 * dist_from)) + cds_sequence + ''.join('N' * dist_to)

                # filter sequence and count vector by AA islands
                count_vector_islands = list()
                cds_sequence_islands = list()

                codon_ind = list()
                for ci in np.array(range(0, len(cds_sequence), 3)):
                    codon = cds_sequence[ci: ci + 3]
                    if codon in aa_codons:
                        codon_ind.append(ci)

                islands = list()
                ii = 0
                while ii < len(codon_ind):
                    if ii == 0:
                        check_left = True
                    else:
                        check_left = codon_ind[ii - 1] - codon_ind[ii] <= dist_from
                    if ii == len(codon_ind) - 1:
                        check_right = True
                    else:
                        check_right = codon_ind[ii] - codon_ind[ii + 1] <= dist_from
                    if check_left and check_right:
                        islands.append(codon_ind[ii])
                        for d in range(dist_from, dist_to):
                            aa_island_count_df.at[aa, d] += count_vector[codon_ind[ii] + d]
                    ii += 1

            self.amino_acid_island_count_df = aa_island_count_df

        return


    # TODO a modification for di-codon counts to me incorporated inseat of get_codon_pauses in the future
    @preconditions(lambda dist_from: isinstance(dist_from, int),
                   lambda dist_from: dist_from < 0,
                   lambda dist_to: isinstance(dist_to, int),
                   lambda dist_to: dist_to >= 0)
    def get_out_of_frame_codon_pauses(self, dist_from=-30, dist_to=3, downsample=True, fstart = 0):
        """
        Counts the meta-number of 5' mapping positions at the given distance from a codon or codon-pair
        Only transcripts with cds of length multiple of 3 are accounted for.
        The only frame in these transcripts is considered.

        :param codon:
        :param dist_from: negative distance from each codon or codon-pair
        :param dist_to: positive distance after each codon or codon-pair
        :param mask_dist: the number of positions to mask in the beginning and end of the gene body
        :param fstart: for -1 frame, fstart = 2; for +1 it is 1
        :return:
        """

        self.logger.info(
            "Counting out of frame codon specific pauses within %d to %d nt distance from the first nucleotide of each codon" %
            (dist_from, dist_to))

        if self.config.args.no_mask:
            mask_dist = 0
            self.logger.info("Transcript boundaries will not be masked")
        else:
            if hasattr(config.args, "codon_mask_size"):
                mask_dist = config.args.codon_mask_size
            else:
                mask_dist = self.MASK_DIST
            self.logger.info("Transcript boundaries will be masked by %d nucleotides" % mask_dist)

        codon_count_df = pd.DataFrame(data=0, index=Codons.CODON_TABLE.keys(),
                                      columns=range(dist_from, dist_to))

        counter = 1

        transcript_assembly = self.annotation.get_transcript_assembly(
            span_size=0)  # don't take more than the gene body (different from previous versions)
        transcript_count = len(transcript_assembly)
        for t in range(transcript_count):
            transcript = transcript_assembly[t]
            if np.floor(transcript_count / 1000) > 0 and counter % 1000 == 0:
                self.logger.info("\r>>Transcript count: %d (%d%s)\t" % (
                    counter, floor(100 * (counter - 1) / transcript_count), '%',), )
            counter += 1

            count_vector = self.get_count_vector(transcript, span_size=0,
                                                 region=FivePSeqCounts.FULL_LENGTH,
                                                 downsample=downsample)

            cds_sequence = self.get_cds_sequence_safe(transcript, 0)

            if sum(count_vector) == 0:
                continue

            if len(cds_sequence) != len(count_vector):
                self.logger.warning("Transcript num %d: cds sequence length %d not equal to count vector length %d"
                                    % (counter, len(cds_sequence), len(count_vector)))
                continue

            # mask the first and last 20 counts to avoid initiation affecting codon-specific counts
            if mask_dist >= 3:
                count_vector[0:mask_dist] = [0] * mask_dist
                cds_sequence = cds_sequence[0:len(cds_sequence) - mask_dist] + ''.join(
                    'N' * (mask_dist - 3)) + cds_sequence[len(cds_sequence) - 3:len(cds_sequence)]

            count_vector = [0] * (-1 * dist_from) + count_vector + [0] * dist_to
            cds_sequence = ''.join('N' * (-1 * dist_from)) + cds_sequence + ''.join('N' * dist_to)

            # identify 3nt bins with non-zero counts
            ind = np.array(range(fstart, len(count_vector), 3))
            hits = [sum(count_vector[i:i + 3]) > 0 for i in ind]
            non_empty_ind = ind[hits]

            # loop through non-empty triplets only
            for i in non_empty_ind:
                # loop through all codons dist_from nucleotides downstream and dist_to nucleotides upstream
                j_range = list(np.arange(i, i - dist_to, -3))[::-1] + list(np.arange(i + 3, i + 3 - dist_from, 3))
                for j in j_range:
                    if j < 0:
                        continue
                    if j + 3 > len(cds_sequence):
                        break
                    codonA = cds_sequence[j: j + 3].upper()

                    if (len(codonA) == 3) & (codonA in Codons.CODON_TABLE.keys()):
                        for p in range(0, 3):
                            d = i - j + p
                            try:
                                codon_count_df.at[codonA, d] += count_vector[i + p]
                            except Exception as e:
                                self.logger.warn("Index out of range: i: %d, j: %d, p: %d, d: %d. %s"
                                                 % (i, j, p, d, str(e)))

        # rename codon_count_df indices by adding amino acid names
        new_index = [Codons.CODON_TABLE.get(codon) + '_' + codon for codon in codon_count_df.index]
        codon_count_df.index = new_index
        codon_count_df = codon_count_df

        return codon_count_df

    def codon_to_amino_acid_count_df(self, codon_count_df):

        amino_acid_count_df = pd.DataFrame(data=0, index=Codons.AMINO_ACID_TABLE.keys(),
                                           columns=codon_count_df.columns)

        for codon in codon_count_df.index:
            aa = Codons.CODON_TABLE.get(codon)
            amino_acid_count_df.loc[aa, :] += codon_count_df.loc[codon, :]

        return amino_acid_count_df

    def get_tripeptide_pos(self):

        if hasattr(config.args, "tripeptide_pos"):
            pos = config.args.tripeptide_pos
        else:
            pos = self.TRIPEPTIDE_POS

        return pos

    def get_dipeptide_pos(self):

        if hasattr(config.args, "dipeptide_pos"):
            pos = config.args.dipeptide_pos
        else:
            pos = self.DIPEPTIDE_POS

        return pos

    def get_cspan(self):

        if hasattr(config.args, "cspan"):
            pos = config.args.cspan
        else:
            pos = self.CSPAN

        return pos


    def filter_codon_counts(self, codon_count_df, pos, top=50):
        """
        Filter the di/tricodon (or di/tripeptide) counts to exclude low counts (rowSums less than the specified threshold) and
        to include only the top di/tricodons with highest relative counts at the given position

        :param codon_count_df: the codon_df to filter
        :param top: the number of highest relative count tricodons to keep
        :param pos: the position to filter the top counts
        :return codon_filtered_df: the filtered count dataframe
        """

        self.logger.info("Sorting and selecting top %d peptides/codons at position %d from the A site" %
                         (top, pos))

        codon_filtered_df = codon_count_df[codon_count_df.sum(1) >= self.COUNT_THRESHOLD]
        pos_rel_counts = codon_filtered_df[pos] / codon_filtered_df.sum(1)
        codon_filtered_df = codon_filtered_df.iloc[
            sorted(range(len(pos_rel_counts)), reverse=True, key=lambda k: pos_rel_counts[k])[0:top]]
        return codon_filtered_df

    def get_amino_acid_stats(self):
        if self.amino_acid_stats_df is None:
            self.amino_acid_stats_df = self.compute_codon_stats_amino_acid()

        return self.amino_acid_stats_df

    def get_codon_stats(self):
        if self.codon_stats_df is None:
            self.codon_stats_df = self.compute_codon_stats_codon()
        return self.codon_stats_df

    def compute_codon_genome_usage(self):
        self.codon_genome_usage_df = pd.DataFrame(data=0, index=Codons.CODON_TABLE.keys(),
                                                  columns=['abs', 'fraction'])
        self.amino_acid_genome_usage_df = pd.DataFrame(data=0, index=Codons.AMINO_ACID_TABLE.keys(),
                                                       columns=['abs', 'fraction'])

    def compute_codon_stats_amino_acid(self):
        return self.compute_codon_stats(self.get_amino_acid_pauses(), self.amino_acid_genome_usage_df)

    def compute_codon_stats_codon(self):
        return self.compute_codon_stats(self.get_codon_pauses(), self.codon_genome_usage_df)

    def compute_codon_stats(self, codon_counts, codon_genome_usage, until=-3):
        """
        Counts usage and frame protection stats for each codon/amino-acid.

        The following dataframe will be generated based on codon counts table:

        codon/aminoacid FPI Frame   peak(pos)   peak(scale) usage(sum of counts)    genome_presence

        :return: dataframe
        """

        self.logger.info("Counting codon usage statistics")

        try:
            stop_ind = codon_counts.keys().to_list().index(until)
            codon_counts = codon_counts.iloc[:, 0:stop_ind]
            f2 = sum([codon_counts.iloc[:, i] for i in reversed(range(stop_ind - 1, -1, -3))])
            f1 = sum([codon_counts.iloc[:, i] for i in reversed(range(stop_ind - 2, -1, -3))])
            f0 = sum([codon_counts.iloc[:, i] for i in reversed(range(stop_ind - 3, -1, -3))])

            codon_stats = pd.DataFrame(list(zip(f0, f1, f2)), columns=['F0', 'F1', 'F2'])

            codon_stats['FPI'] = np.zeros(len(codon_stats))
            codon_stats['F'] = np.zeros(len(codon_stats))
            codon_stats['F_perc'] = np.zeros(len(codon_stats))

            for i in range(len(codon_stats)):
                fpi, fmax, fperc = CountManager.fpi_stats_from_frame_counts(codon_stats.iloc[i, :])
                codon_stats.loc[i, 'FPI'] = fpi
                codon_stats.loc[i, 'F'] = fmax
                codon_stats.loc[i, 'F_perc'] = fperc

            codon_stats['peak_pos'] = [np.argmax(codon_counts.iloc[i, :]) for i in range(len(codon_stats))]
            codon_stats['peak_scale'] = np.zeros(len(codon_stats))

            for i in range(len(codon_stats)):
                for i in range(len(codon_stats)):
                    counts = list(codon_counts.iloc[i, :])
                    if sum(counts) > 0:
                        frame = int(codon_stats.loc[i, 'F'])
                        frame_inds = [j for j in reversed(range(len(counts) - 3 + frame, -1, -3))]
                        frame_counts = [counts[j] for j in frame_inds]
                        codon_stats.loc[i, 'peak_scale'] = len(frame_counts) * max(frame_counts) / sum(frame_counts)
                        codon_stats.loc[i, 'peak_pos'] = codon_counts.columns[frame_inds[np.argmax(frame_counts)]]

            codon_stats['usage'] = list(sum([codon_counts.iloc[:, i] for i in range(0, stop_ind)]))
            codon_stats['genome_usage_abs'] = list(codon_genome_usage.loc[:, 'abs'])
            codon_stats['genome_usage_fraction'] = list(codon_genome_usage.loc[:, 'fraction'])
            usage_norm = codon_stats['usage'] / codon_stats['genome_usage_fraction']
            usage_norm /= sum(usage_norm)
            codon_stats['usage_normalized'] = usage_norm

            codon_stats.index = codon_counts.index

            return codon_stats

        except:
            self.logger.warning("Could not compute codon stats. Codon counts dataframe did not have column %d." % until)
            return None
        # exclude the counts downstream from -3

    @preconditions(lambda loci_file: str)
    def get_pauses_from_loci(self, loci_file, read_locations=READ_LOCATIONS_ALL):
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
            "Counting pauses in %s region from loci given in file %s" % (read_locations, loci_file))

        loci = pd.read_csv(loci_file, sep="\t", index_col=None)
        self.loci_overlaps = []
        # the results will be kept in a dictionary:
        #   key - distance from any locus
        #   value - number of mapping positions at key distance from any locus
        loci_pauses_dict = {}

        span_size = self.annotation.span_size

        counter = 0
        loci_row = 0

        done = False
        move_transcript = True
        move_locus = False
        tg = self.annotation.get_transcript_assembly(span_size)
        transcript = None

        while True:
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

                    if loci.loc[loci_row, "str"] == "+":
                        locus_pos = loci.loc[loci_row, "start"]
                    else:
                        locus_pos = loci.loc[loci_row, "end"]

                    # locus is upstream of transcript -> move locus
                    if transcript.cds_genome_start - span_size > locus_pos:
                        move_locus = True
                        continue
                    # transcript is upstream of locus -> move transcript
                    elif transcript.cds_genome_end + span_size < locus_pos:
                        move_transcript = True
                        continue

                    elif str(transcript.strand) != str(loci.loc[loci_row, "str"]):
                        move_locus = True
                        continue

                    else:
                        count_vector = self.get_count_vector(transcript, span_size, FivePSeqCounts.FULL_LENGTH,
                                                             downsample=True)

                        transcript_genome_start = transcript.cds_genome_start - span_size
                        transcript_genome_end = transcript.cds_genome_end + span_size

                        if len(count_vector) != transcript_genome_end - transcript_genome_start:
                            move_transcript = True
                            continue

                        if transcript.strand == "+":
                            locus_ind = locus_pos - transcript_genome_start
                        else:
                            locus_ind = transcript_genome_end - locus_pos

                        if read_locations == self.READ_LOCATIONS_ALL:
                            ind = np.array(range(len(count_vector) - 2 * span_size, len(count_vector)))
                        elif read_locations == self.READ_LOCATIONS_5UTR:
                            ind = np.array(range(0, span_size))
                        elif read_locations == self.READ_LOCATIONS_3UTR:
                            ind = np.array(range(len(count_vector) - span_size, len(count_vector)))
                        elif read_locations == self.READ_LOCATIONS_CDS:
                            ind = np.array(range(len(count_vector) - 2 * span_size, len(count_vector) - span_size))
                        else:
                            ind = np.array(range(0, len(count_vector)))

                        hits = [count_vector[i] > 0 for i in ind]
                        non_empty_ind = ind[hits]

                        for i in non_empty_ind:
                            distance = i - locus_ind

                            if distance < 2 * span_size and distance >= -2 * span_size:
                                if distance in loci_pauses_dict.keys():
                                    loci_pauses_dict[distance] += count_vector[i]
                                else:
                                    loci_pauses_dict.update({distance: count_vector[i]})

                                overlap = [FivePSeqOut.get_transcript_attr(transcript, "ID"),
                                           FivePSeqOut.get_transcript_attr(transcript, "Name"),
                                           transcript.chrom, transcript.strand,
                                           transcript.cds_genome_start,
                                           transcript.cds_genome_end,
                                           loci.loc[loci_row, "symbol"],
                                           loci.loc[loci_row, "chr"],
                                           loci.loc[loci_row, "str"],
                                           loci.loc[loci_row, "start"],
                                           loci.loc[loci_row, "end"],
                                           i, distance, count_vector[i]]
                                self.loci_overlaps.append(overlap)

                        move_locus = True

                elif str(transcript.chrom) > str(loci.loc[loci_row, "chr"]):
                    move_locus = True
                    continue

                else:
                    move_transcript = True
                    continue
            else:
                break

        # turn the dictionary into a metacount vector, with indices from -1*maxdistance to maxdistance
        self.logger.debug("Merging the dictionary into metacounts")
        maxdist = 2 * span_size
        metacount_vector = [0] * 2 * maxdist
        for i in range(-1 * maxdist, maxdist):
            if i in loci_pauses_dict.keys():
                metacount_vector[maxdist + i] = loci_pauses_dict[i]
        metacount_series = pd.Series(data=metacount_vector, index=np.arange(-1 * maxdist, maxdist))

        return metacount_series

    def get_loci_overlaps_df(self):
        """
        Returns the overlaps of transcripts with given loci in the form of a data-frame with column names.

        :return:
        """
        colnames = ["ID", "Name", "chr", "str", "genome_start", "genome_end", "RBP", "loc_chr", "loc_str", "loc_start",
                    "loc_end",
                    "i", "dist", "count"]

        outliers_df = pd.DataFrame(self.loci_overlaps, index=None, columns=colnames)
        return outliers_df

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

    def filter_pattern(self, transcript, count_vector, span_size):
        """
        Filters out counts proceeding the pattern given in the self.fpat argument.
        Returns the new count vector.

        :param transcript: transcript instance
        :param count_vector: [int] array for counts
        :param span_size: int the span size of the transcript and the count vector (should co-inside)
        :return: [int] array for filtered out counts
        """

        if self.fpat is None:
            pass
        else:
            seq = self.get_sequence(transcript, span_size, span_size)
            if len(seq) != len(count_vector):
                self.logger.warning("Could not filter transcript %s. Reason: different sequence "
                                 "and count vector lengths." % transcript.get_name())
            else:
                ind = np.arange(len(count_vector))
                hits = [count_vector[i] > 0 for i in ind]
                non_empty_ind = ind[hits]
                for i in non_empty_ind:
                    if i - len(self.fpat) < 0:
                        pass
                    elif seq[i - len(self.fpat):i] == self.fpat:
                        count_vector[i] = 0
        return count_vector

class FivePSeqCountsContainer:
    """
    A wraper for the following data structures:
        count_vector_list_start = None
        count_vector_list_term = None
        count_vector_list_full_length = None
        meta_count_series_start = None
        meta_count_series_term = None
        frame_counts_df_start = None
        frame_counts_df_term = None

    """
    count_vector_list_start = None
    count_vector_list_term = None
    count_vector_list_full_length = None
    meta_count_series_start = None
    meta_count_series_term = None
    frame_counts_df_start = None
    frame_counts_df_term = None

    def __init__(self, count_vector_list_start,
                 count_vector_list_term,
                 count_vector_list_full_length,
                 meta_count_series_start,
                 meta_count_series_term,
                 frame_counts_df_start,
                 frame_counts_df_term):
        self.count_vector_list_start = count_vector_list_start
        self.count_vector_list_term = count_vector_list_term
        self.count_vector_list_full_length = count_vector_list_full_length
        self.meta_count_series_term = meta_count_series_term
        self.meta_count_series_start = meta_count_series_start
        self.frame_counts_df_start = frame_counts_df_start
        self.frame_counts_df_term = frame_counts_df_term


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
        max_len = 0
        for i in range(len(count_vector_list)):
            if len(count_vector_list) > max_len:
                max_len = len(count_vector_list[i])

        for i in range(len(count_vector_list)):
            if len(count_vector_list[i]) < max_len:
                short_vec = count_vector_list[i]
                long_vec = [0] * max_len
                long_vec[0:len(short_vec)] = short_vec
                count_vector_list[i] = long_vec

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
            logger = logging.getLogger(config.FIVEPSEQ_LOGGER)
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
            logger = logging.getLogger(config.FIVEPSEQ_LOGGER)
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

        elif region == FivePSeqCounts.MID:
            d = np.arange(0, len(count_vector))

        else:
            error_message = "Invalid region %s specified: should be either %s or %s" \
                            % (region, FivePSeqCounts.START, FivePSeqCounts.TERM)
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
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

        logging.getLogger(config.FIVEPSEQ_LOGGER).debug(
            "Retrieving count-sums per frame relative to %s ..." % region)

        # Create an empty dataframe
        n = len(count_vector_list)
        frame_counts_df = pd.DataFrame({'F0': [0] * n,
                                        'F1': [0] * n,
                                        'F2': [0] * n})
        for t_ind in range(0, n):
            # Print status update to console
            if t_ind % 10000 == 0:
                logging.getLogger(config.FIVEPSEQ_LOGGER).info("\r>>Transcript count: %d (%d%s)\t" % (
                    t_ind, np.floor(100 * (t_ind - 1) / n), '%'))

            # extract frame count vectors from count vectors
            count_vector = count_vector_list[t_ind]
            if sum(count_vector) == 0:
                for f in range(0, 3):
                    frame_counts_df.iloc[t_ind, f] = 0
            else:
                frame_counts = CountManager.extract_frame_count_vectors(count_vector, span_size, region)

                # sum-up counts in each frame and add to the dataframe
                for f in range(0, 3):
                    frame_counts_df.iloc[t_ind, f] = sum(frame_counts[f])

        return frame_counts_df

    @staticmethod
    @preconditions(lambda file_path: isinstance(file_path, str))
    def read_index_as_list(file_path):
        """
        Reads a new line separated list of integers from a file to a list of indices.

        :param file_path: str: full path to the file
        :return: [int]: list of indices
        :exception: raise IOError if file does not exist
        """
        if not os.path.exists(file_path):
            error_message = "Problem reading counts: the file %s does not exist" % file_path
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_LOGGER).debug("Reading count file %s" % file_path)

        indices = list(pd.read_csv(file_path, header=None).iloc[:, 0])
        return indices

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
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_LOGGER).debug("Reading count file %s" % file_path)
        df = pd.read_csv(file_path, header=None, sep="|")
        count_vector_list = [[]] * len(df)
        for i in range(0, len(df)):
            count_vector_list[i] = list(map(int, df.iloc[i, 0].split("\t")))
        return count_vector_list

    @staticmethod
    @preconditions(lambda count_vec_list: isinstance(count_vec_list, list))
    def extract_mid_counts(count_vec_list, size=600, tail=100):
        """
        Extract the middle parts of each vector in the given list
        and returns the list of extracted same-length vectors.

        The frame is preserved and adjusted by the tail length.

        :param count_vec_list: list of full size vectors
        :param size: the size of the middle part to extract
        :param tail: the distance spanning from CDS start and end in the given list of vectors
        :return: a list of frame-adjusted length sized vectors from the middle part of each vector
        """

        mid_vec_list = [[]] * len(count_vec_list)

        size = size - size % 3

        for i in range(len(count_vec_list)):
            vec = count_vec_list[i]
            # determine start, putting 0 at the first nucleotide of the start codon
            cds_vec = vec[0 + tail: len(vec) - tail - 3: 1]  # TODO check this

            if (len(cds_vec) < size):
                mid_vec = [0] * int(floor((size - len(cds_vec)) / 6) * 3) + \
                          cds_vec + \
                          [0] * int(np.ceil((size - len(cds_vec)) / 6) * 3)

            else:
                num_f0_bases = len(cds_vec[0: len(cds_vec): 3])
                mid_start = int(np.ceil(num_f0_bases - size / 3) / 2) * 3
                mid_vec = cds_vec[mid_start:mid_start + size]

            mid_vec_list[i] = mid_vec

        return mid_vec_list

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
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_LOGGER).debug("Reading meta count file %s" % file)
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
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_LOGGER).debug("Reading frame counts file %s" % file)

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
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_LOGGER).debug("Reading amino acids pauses file %s" % file)

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
            logging.getLogger(config.FIVEPSEQ_LOGGER).debug(
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
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(error_message)
            raise IOError(error_message)
        logging.getLogger(config.FIVEPSEQ_LOGGER).debug("Reading count file %s" % file_path)
        if os.stat(file_path).st_size == 0:
            counts = []
        else:
            logging.getLogger("Reading in count distribution (this may last a few minutes for large libraries)")
            counts = pd.read_csv(file_path, header=None)
            counts = list(map(int, counts.iloc[:, 0]))

        return counts

    @staticmethod
    @preconditions(lambda file_path: isinstance(file_path, str))
    def read_outlier_lower(file_path):
        """
        Reads and returns the lower value of read numbers to be considered as outliers for downsampling.

        :param file_path: a file containing a single float number
        :return: float
        """

        outlier_lower = float(pd.read_csv(file_path, header=None))
        return outlier_lower

    @staticmethod
    @preconditions(lambda file_path: isinstance(file_path, str))
    def read_count_dict(file_path):
        """
        Reads a tab-delimited file and returns a dictionary of count frequencies.

        :param file_path: a file containing a dictionary of count frequencies
        :return: float
        """
        count_freq_dict = {}
        dict_mat = pd.read_csv(file_path, header=None, delimiter="\t", index_col=0)
        for i in range(len(dict_mat)):
            count_freq_dict[dict_mat.index[i]] = dict_mat.iloc[i, 0]

        return collections.OrderedDict(sorted(count_freq_dict.items()))

    @staticmethod
    def filter_fivepseqCountsContainer(fivepseqcountsContainer, transcript_indices, span_size=100):
        """
        Gets a fivepseq_counts instance with non-empty count items and filters each by provided transcript indices
            count_vector_list_start = None
            count_vector_list_term = None
            count_vector_list_full_length = None
            meta_count_series_start = None
            meta_count_series_term = None
            frame_counts_df_start = None
            frame_counts_df_term = None

        :param fivepseqcountsContainer:
        :param transcript_ind:
        :return:
        """
        if fivepseqcountsContainer.count_vector_list_full_length is not None:
            count_vector_list_full_length = [fivepseqcountsContainer.count_vector_list_full_length[i] for i
                                             in transcript_indices]
        else:
            count_vector_list_full_length = None

        count_vector_list_term = [fivepseqcountsContainer.count_vector_list_term[i] for i in
                                  transcript_indices]

        count_vector_list_start = [fivepseqcountsContainer.count_vector_list_start[i] for i in
                                   transcript_indices]

        meta_count_series_term = CountManager.count_vector_to_df(
            CountManager.compute_meta_counts(count_vector_list_term),
            FivePSeqCounts.TERM, tail=span_size)

        meta_count_series_start = CountManager.count_vector_to_df(
            CountManager.compute_meta_counts(count_vector_list_start),
            FivePSeqCounts.START, tail=span_size)

        frame_counts_df_term = fivepseqcountsContainer.frame_counts_df_term.iloc[transcript_indices,]

        frame_counts_df_start = fivepseqcountsContainer.frame_counts_df_start.iloc[
            transcript_indices,]

        fivepseq_counts_filtered = FivePSeqCountsContainer(count_vector_list_start, count_vector_list_term,
                                                           count_vector_list_full_length,
                                                           meta_count_series_start, meta_count_series_term,
                                                           frame_counts_df_start, frame_counts_df_term)

        return fivepseq_counts_filtered

    @staticmethod
    def combine_count_series(count_series_dict, lib_size_dict=None, scale=False):
        """
        Combines counts in the series dictionary and returns a single count series.
        If lib_size_dict is not None, than the counts are first weighted based on the library size and then combined.
        Weighting is done in a way to give higher weight to samples with larger library sizes.

        :param count_series_dict:
        :param lib_size_dict:
        :return:
        """
        count_series_combined = None
        start = True
        for key in count_series_dict.keys():
            count_series = count_series_dict[key].copy()
            if lib_size_dict is not None:
                if scale:
                    count_series.C /= (float(lib_size_dict[key]) / (10 ** 6)) * len(lib_size_dict)
                else:
                    count_series.C *= float(lib_size_dict[key]) / sum(lib_size_dict.values())

            if start:
                count_series_combined = count_series.copy()
                start = False
            else:
                count_series_combined.C += count_series.C

        return count_series_combined

    @staticmethod
    def combine_frame_counts(frame_count_dict, lib_size_dict=None):
        """
        Combines counts in the dataframe dictionary and returns a single dataframe.
        If lib_size_dict is not None, than the counts are first weighted based on the library size and then combined.
        Weighting is done in a way to give higher weight to samples with larger library sizes.

        :param count_series_dict:
        :param lib_size_dict:
        :return:
        """
        frame_count_combined = None
        start = True
        for key in frame_count_dict.keys():
            count_df = frame_count_dict[key]
            if lib_size_dict is not None:
                count_df.loc[:, ('F0', 'F1', 'F2')] *= float(lib_size_dict[key]) / sum(lib_size_dict.values())

            if start:
                frame_count_combined = count_df.copy()
                start = False
            else:
                frame_count_combined.loc[:, ('F0', 'F1', 'F2')] += count_df.loc[:, ('F0', 'F1', 'F2')]

        return frame_count_combined

    @staticmethod
    def combine_amino_acid_dfs(amino_acid_df_dict, lib_size_dict=None):
        """
        Combines counts in the dataframe dictionary and returns a single dataframe.
        If lib_size_dict is not None, than the counts are first weighted based on the library size and then combined.
        Weighting is done in a way to give higher weight to samples with larger library sizes.

        :param count_series_dict:
        :param lib_size_dict:
        :return:
        """
        amino_acid_df_combined = None
        start = True
        for key in amino_acid_df_dict.keys():
            count_df = amino_acid_df_dict[key].copy(deep=True)
            if lib_size_dict is not None:
                count_df *= float(lib_size_dict[key]) / sum(lib_size_dict.values())

            if start:
                amino_acid_df_combined = count_df.copy()
                start = False
            else:
                amino_acid_df_combined += count_df

        return amino_acid_df_combined

    @staticmethod
    def fpi_stats_from_frame_counts(frame_counts):
        """
        Takes as input a vector named [F0, F1, F2] and returns:
        (fpi, fmax, f_perc)
        fpi = frame protection index of the maximum frame
        fmax = the maximum frame
        f_perc = the fraction of counts in the maximum frame

        :param frame_counts:
        :return:
        """
        f_counts = (frame_counts['F0'], frame_counts['F1'], frame_counts['F2'])
        fmax = np.argmax(f_counts)
        nom = f_counts[fmax]
        if nom == 0:
            fpi = None
            f_perc = None
        else:
            denom = (sum(f_counts) - nom) / 2.
            if denom == 0:
                fpi = np.log2(float(nom) / 0.5)
            else:
                fpi = np.log2(float(nom / denom))
            f_perc = 100 * (float(f_counts[fmax]) / sum(f_counts))

        return fpi, fmax, f_perc
