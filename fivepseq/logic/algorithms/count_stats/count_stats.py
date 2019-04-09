import logging
import os

import numpy as np
import pandas as pd
from scipy.stats import stats

from fivepseq import config
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts, CountManager
from fivepseq.util.writers import FivePSeqOut


class CountStats:
    fivepseq_counts = None
    fivepseq_out = None
    config = None
    logger = None

    F0 = "F0"
    F1 = "F1"
    F2 = "F2"

    FRAME_COUNT = "f_count"
    FRAME_PERC = "f_perc"
    FPI = "fpi"
    PVAL_PAIR = "p_val_pair"
    PVAL_PAIR_MAX = "p_val_pair_max"
    PVAL_FPI = "p_val_fpi"

    # total stat strings
    TOTAL_NUM_READS = "NumOfReads"
    TOTAL_NUM_POSITIONS = "NumOfMapPositions"
    TOTAL_NUM_TRANSCRIPTS = "NumOfTranscripts"
    TOTAL_NUM_MAP_TRANSCRIPTS = "NumOfMapTranscripts"
    MIN_NUM_READS_PER_TRANSCRIPT = "MinNumOfReadsPerTranscript"
    MAX_NUM_READS_PER_TRANSCRIPT = "MaxNumOfReadsPerTranscript"
    MEDIAN_NUM_READS_PER_TRANSCRIPT = "MedianNumOfReadsPerTranscript"

    data_summary_series = None
    # total num
    total_num_reads = 0
    total_num_positions = 0
    total_num_map_transcripts = 0
    min_num_reads_per_transcript = 0
    max_num_reads_per_transcript = 0
    median_num_reads_per_transcript = 0

    # dfs
    frame_stats_df = None
    fft_stats_start = None
    fft_stats_term = None

    def __init__(self, fivepseq_counts, fivepseq_out, config):
        self.fivepseq_counts = fivepseq_counts
        self.fivepseq_out = fivepseq_out
        self.config = config
        self.logger = logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER)

    def count_stats(self):
        """
        Returns the statistics on frame-related periodicity from a list of count vectors.

        :return:
        """
        self.summarize_transcript_stats()

        if self.data_summary_series[self.TOTAL_NUM_READS] == 0:
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).warning("Zero reads in coding regions: frame and fft stats will not be calculated")
        else:
            self.compute_frame_preference_stats()

            self.compute_fft_stats()

    def summarize_transcript_stats(self):
        self.logger.info("Summarizing transcript counts")
        transcript_descriptors_f = self.fivepseq_out.get_file_path(FivePSeqOut.TRANSCRIPT_DESCRIPTORS_FILE)

        if not os.path.exists(transcript_descriptors_f):
            self.logger.warning("Transcript descriptors file %s not found. Data summary file will not be generated"
                                % transcript_descriptors_f)
        else:
            if os.path.exists(self.fivepseq_out.get_file_path(FivePSeqOut.DATA_SUMMARY_FILE)):
                self.logger.info("Using existing file %s" % FivePSeqOut.DATA_SUMMARY_FILE)
                self.data_summary_series = pd.read_csv(self.fivepseq_out.get_file_path(FivePSeqOut.DATA_SUMMARY_FILE), sep = "\t",  header = None, index_col = 0).iloc[:,0]
            else:
                self.data_summary_series = pd.Series(name=(self.TOTAL_NUM_READS,
                                                           self.TOTAL_NUM_POSITIONS,
                                                           self.TOTAL_NUM_TRANSCRIPTS,
                                                           self.TOTAL_NUM_MAP_TRANSCRIPTS,
                                                           self.MIN_NUM_READS_PER_TRANSCRIPT,
                                                           self.MAX_NUM_READS_PER_TRANSCRIPT,
                                                           self.MEDIAN_NUM_READS_PER_TRANSCRIPT))
                self.logger.debug("Reading file %s" % transcript_descriptors_f)
                transcript_descriptors = pd.read_csv(transcript_descriptors_f, sep="\t", index_col=0)

                self.data_summary_series[self.TOTAL_NUM_READS] = int(sum(
                    transcript_descriptors.loc[:, FivePSeqCounts.NUMBER_READS]))
                self.data_summary_series[self.TOTAL_NUM_POSITIONS] = int(sum(
                    transcript_descriptors.loc[:, FivePSeqCounts.NUMBER_POSITIONS]))
                self.data_summary_series[self.TOTAL_NUM_TRANSCRIPTS] = int(transcript_descriptors.shape[0])
                self.data_summary_series[self.TOTAL_NUM_MAP_TRANSCRIPTS] = int(
                    np.count_nonzero(transcript_descriptors.loc[:, FivePSeqCounts.NUMBER_READS]))
                self.data_summary_series[self.MIN_NUM_READS_PER_TRANSCRIPT] = int(np.min(
                    transcript_descriptors.loc[:, FivePSeqCounts.NUMBER_READS]))
                self.data_summary_series[self.MAX_NUM_READS_PER_TRANSCRIPT] = int(np.max(
                    transcript_descriptors.loc[:, FivePSeqCounts.NUMBER_READS]))
                self.data_summary_series[self.MEDIAN_NUM_READS_PER_TRANSCRIPT] = int(np.median(
                    transcript_descriptors.loc[:, FivePSeqCounts.NUMBER_READS]))

                self.logger.info("Writing summary data")
                self.fivepseq_out.write_series_to_file(self.data_summary_series, FivePSeqOut.DATA_SUMMARY_FILE)

    def compute_fft_stats(self):
        if config.args.conflicts == config.ADD_FILES and \
                (os.path.exists(self.fivepseq_out.get_file_path(FivePSeqOut.FFT_SIGNALS_TERM)) &
                os.path.exists(self.fivepseq_out.get_file_path(FivePSeqOut.FFT_SIGNALS_START)) &
                os.path.exists(self.fivepseq_out.get_file_path(FivePSeqOut.FFT_STATS_DF_FILE))):
            self.logger.info("Skipping FFT statistics calculation: files already exist")
        else:
            self.logger.info("Computing FFT statistics")
            count_vector_list = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH)

            span_size = self.fivepseq_counts.annotation.span_size

            lengths = [0] * len(count_vector_list)
            for i in range(len(count_vector_list)):
                count_vector = count_vector_list[i][span_size:len(count_vector_list[i]) - span_size]
                lengths[i] = len(count_vector)

            # align start
            size = int(stats.scoreatpercentile(lengths, per=25))
            num = 3 * len(count_vector_list) // 4
            start_array = np.zeros((num, size))
            term_array = np.zeros((num, size))

            ind = 0
            for i in range(len(count_vector_list)):
                count_vector = count_vector_list[i][span_size:len(count_vector_list[i]) - span_size]
                if (len(count_vector)) > size:
                    start_array[ind, :] = count_vector[0:size]
                    term_array[ind, :] = count_vector[len(count_vector) - size:len(count_vector)]
                    ind += 1
            start_mean = start_array.mean(axis=0)
            term_mean = term_array.mean(axis=0)

            self.fft_stats_start = self.fft_stats_on_vector(start_mean, 5)
            self.fft_stats_term = self.fft_stats_on_vector(term_mean, 5)
            self.fft_stats_df = pd.DataFrame(
                data=[self.fft_stats_start[1], self.fft_stats_start[2], self.fft_stats_start[3],
                      self.fft_stats_term[1], self.fft_stats_term[2], self.fft_stats_term[3]],
                index=["START_periods", "START_signals", "START_scales",
                       "TERM_periods", "TERM_signals", "TERM_scales"]
            )
            self.logger.info("Writing FFT stats")
            self.fivepseq_out.write_df_to_file(self.fft_stats_df, FivePSeqOut.FFT_STATS_DF_FILE)
            self.fivepseq_out.write_series_to_file(
                self.fft_stats_start[0],
                FivePSeqOut.FFT_SIGNALS_START)
            self.fivepseq_out.write_series_to_file(
                self.fft_stats_term[0],
                FivePSeqOut.FFT_SIGNALS_TERM)

    @staticmethod
    def fft_stats_on_vector(count_vector, n=1):
        fft_abs = np.abs(np.fft.fft(count_vector))
        fft_abs = fft_abs[1:(len(fft_abs) // 2 + 1)]

        ind = np.flip(np.argsort(fft_abs)[-1 * n:])
        signals = [fft_abs[i] for i in ind]
        scales = [signal / np.mean(fft_abs) for signal in signals]
        periods = [0] * n
        ii = 0
        for i in ind:
            period = float(len(count_vector) + 1) / (i + 1)
            if period > len(count_vector) // 2:
                period = 1
            periods[ii] = period
            ii += 1

        d = [float(len(count_vector) + 1) / (i + 1) for i in range(len(fft_abs))]
        fft_signal_series = pd.Series(data=fft_abs, index=d)

        return fft_signal_series, periods, signals, scales

    def compute_frame_preference_stats(self):
        """
        percentage_per_frame: F0, F1, F2
        frame_protection_index: F0, F1, F2
        p-values: (F0-F1, F1-F2, F2-F0)

        :param frame_counts_df:
        :return: frame_stats_df
        """
        if config.args.conflicts == config.ADD_FILES and os.path.exists(self.fivepseq_out.get_file_path(FivePSeqOut.FRAME_STATS_DF_FILE)):
            self.logger.info("Skipping frame stats calculation: file %s exists" % FivePSeqOut.FRAME_STATS_DF_FILE)
        else:
            frame_counts_df = self.fivepseq_counts.get_frame_counts_df(FivePSeqCounts.START)
            self.logger.info("Computing frame preference statistics")
            total_counts = sum(frame_counts_df.iloc[:, 0]) \
                           + sum(frame_counts_df.iloc[:, 1]) \
                           + sum(frame_counts_df.iloc[:, 2])

            frame_names = [self.F0, self.F1, self.F2]
            frame_stats_df = pd.DataFrame(index=[self.FRAME_COUNT, self.FRAME_PERC, self.FPI,
                                                 self.PVAL_PAIR, self.PVAL_FPI,
                                                 self.PVAL_PAIR_MAX],
                                          columns=frame_names)

            if total_counts == 0:
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).warning("Total counts on frames equal to 0")
            else:

                for i in range(3):
                    frame_stats_df.loc[self.FRAME_COUNT, frame_names[i]] = int(sum(frame_counts_df.iloc[:, i]))

                    frame_stats_df.loc[self.FRAME_PERC, frame_names[i]] = 100 * float(
                        sum(frame_counts_df.iloc[:, i])) / total_counts

                    if (total_counts - sum(frame_counts_df.iloc[:, i])) == 0:
                        frame_stats_df.loc[self.FPI, frame_names[i]] = -1
                    else:
                        frame_stats_df.loc[self.FPI, frame_names[i]] = np.log2(
                            float(sum(frame_counts_df.iloc[:, i])) /
                            ((total_counts - sum(frame_counts_df.iloc[:, i])) / 2))

                    frame_stats_df.loc[self.PVAL_PAIR, frame_names[i]] = stats.ttest_ind \
                        (frame_counts_df.iloc[:, i],
                         frame_counts_df.iloc[:, (i + 1) % 3]
                         ).pvalue

                    frame_stats_df.loc[self.PVAL_FPI, frame_names[i]] = stats.ttest_ind \
                        (frame_counts_df.iloc[:, i],
                         list(frame_counts_df.iloc[:, (i + 1) % 3]) + list(frame_counts_df.iloc[:, (i + 2) % 3])
                         ).pvalue

                for i in range(3):
                    frame_stats_df.loc[self.PVAL_PAIR_MAX, frame_names[i]] = np.max(
                        [frame_stats_df.loc[self.PVAL_PAIR, :][(i - 1) % 3],
                         frame_stats_df.loc[self.PVAL_PAIR, :][i]])

            self.frame_stats_df = frame_stats_df

            self.fivepseq_out.write_df_to_file(self.frame_stats_df, FivePSeqOut.FRAME_STATS_DF_FILE)


    def get_frame_of_preference(self):
        if self.frame_stats_df is None:
            self.compute_frame_preference_stats()

        frame = np.argmax(list(self.frame_stats_df.loc[self.FPI, :]))
        fpi = self.frame_stats_df.iloc[self.FPI, frame]
        pval_fpi = self.frame_stats_df.iloc[self.PVAL_FPI, frame]

        return frame, fpi, pval_fpi

    def get_significant_frame(self):
        if self.frame_stats_df is None:
            self.compute_frame_preference_stats()

        frame = np.argmin(list(self.frame_stats_df.loc[self.PVAL_PAIR_MAX, :]))
        fpi = self.frame_stats_df.iloc[self.FPI, frame]
        pval_pair_max = self.frame_stats_df.iloc[self.PVAL_PAIR_MAX, frame]

        return frame, fpi, pval_pair_max
