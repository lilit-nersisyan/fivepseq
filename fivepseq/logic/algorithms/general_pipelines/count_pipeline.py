import os

import config
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts, CountManager
from logic.algorithms.count_stats.count_stats import CountStats


class CountPipeline:
    fivepseq_counts = None
    fivepseq_out = None

    def __init__(self, fivepseq_counts, fivepseq_out):
        self.fivepseq_counts = fivepseq_counts
        self.fivepseq_out = fivepseq_out

    def run(self, logger):
        #   terminal counts
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_TERM_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.COUNT_TERM_FILE)
        else:
            term_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.TERM)
            self.fivepseq_out.write_vector_list(term_counts, self.fivepseq_out.COUNT_TERM_FILE)

        #   start counts
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_START_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.COUNT_START_FILE)
        else:
            start_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.START)
            self.fivepseq_out.write_vector_list(start_counts, self.fivepseq_out.COUNT_START_FILE)

        #   full length counts
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_FULL_FILE))):
            # NOTE I need this done fast now - so I'll read from file
            self.fivepseq_counts.set_count_vector_list(CountManager.read_counts_as_list(
                self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_FULL_FILE)), FivePSeqCounts.FULL_LENGTH)
            logger.debug("Skipping existing file %s" % self.fivepseq_out.COUNT_FULL_FILE)
        else:
            full_length_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH)
            self.fivepseq_out.write_vector_list(full_length_counts, self.fivepseq_out.COUNT_FULL_FILE)

        #   meta counts
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.META_COUNT_TERM_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.META_COUNT_TERM_FILE)
        else:
            self.fivepseq_out.write_series_to_file(self.fivepseq_counts.get_meta_count_series(FivePSeqCounts.TERM),
                                                   self.fivepseq_out.META_COUNT_TERM_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.META_COUNT_START_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.META_COUNT_START_FILE)
        else:
            self.fivepseq_out.write_series_to_file(self.fivepseq_counts.get_meta_count_series(FivePSeqCounts.START),
                                                   self.fivepseq_out.META_COUNT_START_FILE)
        #   frame counts
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_START_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.FRAME_COUNTS_START_FILE)
        else:
            self.fivepseq_out.write_df_to_file(
                CountManager.extract_count_sums_per_frame_per_transcript(
                    self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH),
                    self.fivepseq_counts.span_size, FivePSeqCounts.START),
                self.fivepseq_out.FRAME_COUNTS_START_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_TERM_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.FRAME_COUNTS_TERM_FILE)
        else:
            self.fivepseq_out.write_df_to_file(
                CountManager.extract_count_sums_per_frame_per_transcript(
                    self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH),
                    self.fivepseq_counts.span_size, FivePSeqCounts.TERM),
                self.fivepseq_out.FRAME_COUNTS_TERM_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.START_CODON_DICT_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.START_CODON_DICT_FILE)
        else:
            self.fivepseq_out.write_dict(self.fivepseq_counts.start_codon_dict,
                                         self.fivepseq_out.START_CODON_DICT_FILE)
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.TERM_CODON_DICT_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.TERM_CODON_DICT_FILE)
        else:
            self.fivepseq_out.write_dict(self.fivepseq_counts.stop_codon_dict,
                                         self.fivepseq_out.TERM_CODON_DICT_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)
        else:
            self.fivepseq_out.write_vector(self.fivepseq_counts.canonical_transcript_index,
                                           self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)
        else:
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.transcript_descriptors,
                                               self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.AMINO_ACID_PAUSES_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.AMINO_ACID_PAUSES_FILE)
        else:
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_amino_acid_pauses(
                self.fivepseq_counts.span_size, 20),
                self.fivepseq_out.AMINO_ACID_PAUSES_FILE)

        #   transcript assembly
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE)
        else:
            self.fivepseq_out.write_transcript_assembly_to_file(
                self.fivepseq_counts.annotation.transcript_assembly,
                self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE)

        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.LOCI_PAUSES_FILE))):
            logger.debug("Skipping existing file %s" % self.fivepseq_out.LOCI_PAUSES_FILE)
        else:
            if self.fivepseq_counts.loci_file is not None:
                self.fivepseq_out.write_vector(
                    self.fivepseq_counts.get_pauses_from_loci(self.fivepseq_counts.loci_file,
                                                              self.fivepseq_counts.span_size,
                                                              20),
                    self.fivepseq_out.LOCI_PAUSES_FILE)

        # frame stats
        frame_stats = CountStats.frame_preference_stats(
            self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH))
