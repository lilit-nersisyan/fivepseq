import numpy as np

class CountStats:

    #strings
    START_CODON = "start"
    STOP_CODON = "stop"
    TRANSCRIPT_LENGTH = "len"
    TRANSCRIPT_3NT = "3nt"
    NUMBER_MAPPINGS = "NumOfReads"
    NUMBER_POSITIONS = "NumOfMapPositions"

    # total stat strings
    TOTAL_NUM_READS = "NumOfReads"
    TOTAL_NUM_POSITIONS = "NumOfMapPositions"
    TOTAL_NUM_MAP_TRANSCRIPTS = "NumOfMapTranscripts"
    MIN_NUM_READS_PER_TRANSCRIPT = "MinNumOfReadsPerTranscript"
    MAX_NUM_READS_PER_TRANSCRIPT = "MaxNumOfReadsPerTranscript"
    MEDIAN_NUM_READS_PER_TRANSCRIPT = "MedianNumOfReadsPerTranscript"

    # total num
    total_num_reads = 0
    total_num_positions = 0
    total_num_map_transcripts = 0
    min_num_reads_per_transcript = 0
    max_num_reads_per_transcript = 0
    median_num_reads_per_transcript = 0

    frame_protection_index = 0
    preferred_frame = 0
    fft_signal = 0
    fft_pval = -1




    @staticmethod
    def frame_preference_stats(count_vector_list):
        """
        Returns the statistics on frame-related periodicity from a list of count vectors.

        :param count_vector_list:
        :return:
        """
        df = pd.read_csv(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_TERM_FILE), sep="\t")
        tt = stats.ttest_ind(df.iloc[:, 1], df.iloc[:, 1])
        for count_vector in count_vector_list:
            np.fft.fft(count_vector, 3)