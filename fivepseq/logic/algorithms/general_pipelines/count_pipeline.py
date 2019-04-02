import logging
import os

from fivepseq.util.writers import FivePSeqOut

import config
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts, CountManager
from logic.algorithms.count_stats.count_stats import CountStats


class CountPipeline:
    fivepseq_counts = None
    fivepseq_out = None

    def __init__(self, fivepseq_counts, fivepseq_out):
        self.fivepseq_counts = fivepseq_counts
        self.fivepseq_out = fivepseq_out

    def run(self):
        logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info("\n\nFivepseq started for bam file %s\n\n"
                                                             % os.path.basename(
            self.fivepseq_counts.alignment.bam_file))

        #   full length counts
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_FULL_FILE)) &
                os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE))):

            # NOTE I need this done fast now - so I'll read from file
            self.fivepseq_counts.set_count_vector_list(CountManager.read_counts_as_list(
                self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_FULL_FILE)), FivePSeqCounts.FULL_LENGTH)

            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                "Skipping existing file %s" % self.fivepseq_out.COUNT_FULL_FILE)
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                "Skipping existing file %s" % self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)
        else:
            self.fivepseq_counts.generate_count_vector_lists()
            full_length_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH)
            self.fivepseq_out.write_vector_list(full_length_counts, self.fivepseq_out.COUNT_FULL_FILE)
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.transcript_descriptors,
                                               self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)

        # frame stats
        count_stats = CountStats(self.fivepseq_counts, self.fivepseq_out)
        count_stats.count_stats()

        if count_stats.data_summary_series[count_stats.TOTAL_NUM_READS] == 0:
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).warning(
                "No reads found in coding regions. Fivepseq will skip the rest of calculations.")

        else:
            #   terminal counts
            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_TERM_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.COUNT_TERM_FILE)
            else:
                term_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.TERM)
                self.fivepseq_out.write_vector_list(term_counts, self.fivepseq_out.COUNT_TERM_FILE)

            #   start counts
            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_START_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.COUNT_START_FILE)
            else:
                start_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.START)
                self.fivepseq_out.write_vector_list(start_counts, self.fivepseq_out.COUNT_START_FILE)

            #   meta counts
            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.META_COUNT_TERM_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.META_COUNT_TERM_FILE)
            else:
                self.fivepseq_out.write_series_to_file(self.fivepseq_counts.get_meta_count_series(FivePSeqCounts.TERM),
                                                       self.fivepseq_out.META_COUNT_TERM_FILE)

            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.META_COUNT_START_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.META_COUNT_START_FILE)
            else:
                self.fivepseq_out.write_series_to_file(self.fivepseq_counts.get_meta_count_series(FivePSeqCounts.START),
                                                       self.fivepseq_out.META_COUNT_START_FILE)
            #   frame counts
            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_START_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.FRAME_COUNTS_START_FILE)
            else:
                frame_count_df_start = self.fivepseq_counts.get_frame_counts_df(FivePSeqCounts.START)
                self.fivepseq_out.write_df_to_file(
                    frame_count_df_start,
                    self.fivepseq_out.FRAME_COUNTS_START_FILE)

            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_TERM_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.FRAME_COUNTS_TERM_FILE)
            else:
                frame_count_df_term = self.fivepseq_counts.get_frame_counts_df(FivePSeqCounts.TERM)
                self.fivepseq_out.write_df_to_file(
                    frame_count_df_term,
                    self.fivepseq_out.FRAME_COUNTS_TERM_FILE)

            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.START_CODON_DICT_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.START_CODON_DICT_FILE)
            else:
                self.fivepseq_out.write_dict(self.fivepseq_counts.start_codon_dict,
                                             self.fivepseq_out.START_CODON_DICT_FILE)
            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.TERM_CODON_DICT_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.TERM_CODON_DICT_FILE)
            else:
                self.fivepseq_out.write_dict(self.fivepseq_counts.stop_codon_dict,
                                             self.fivepseq_out.TERM_CODON_DICT_FILE)

            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)
            else:
                self.fivepseq_out.write_vector(self.fivepseq_counts.canonical_transcript_index,
                                               self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)

            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.AMINO_ACID_PAUSES_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.AMINO_ACID_PAUSES_FILE)
            else:
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_amino_acid_pauses(
                    0, 20),
                    self.fivepseq_out.AMINO_ACID_PAUSES_FILE)

            #   transcript assembly
            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE)
            else:
                self.fivepseq_out.write_transcript_assembly_to_file(
                    self.fivepseq_counts.annotation.transcript_assembly,
                    self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE)

            if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                    os.path.exists(self.fivepseq_out.get_file_path(self.fivepseq_out.LOCI_PAUSES_FILE))):
                logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                    "Skipping existing file %s" % self.fivepseq_out.LOCI_PAUSES_FILE)
            else:
                if self.fivepseq_counts.loci_file is not None:
                    self.fivepseq_out.write_series_to_file(
                        self.fivepseq_counts.get_pauses_from_loci(self.fivepseq_counts.loci_file,
                                                                  self.fivepseq_counts.span_size,
                                                                  20),
                        self.fivepseq_out.LOCI_PAUSES_FILE)

        success = self.fivepseq_out.sanity_check_for_counts()
        if success:
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                "\n\nFivepseq successfully finished for bam file %s\n\n"
                % os.path.basename(self.fivepseq_counts.alignment.bam_file))
        else:
            logging.getLogger(config.FIVEPSEQ_COUNT_LOGGER).info(
                "\n\nFivepseq finished for bam file %s.\n Some files failed to be generated. Check those in %s.\n\n"
                % (os.path.basename(self.fivepseq_counts.alignment.bam_file)),
                self.fivepseq_out.get_file_path(FivePSeqOut.FAILED_COUNT_FILES_LIST))
