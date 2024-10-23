import logging
import os

from fivepseq.util.writers import FivePSeqOut

from fivepseq import config
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts, CountManager
from fivepseq.logic.algorithms.count_stats.count_stats import CountStats
from fivepseq.logic.algorithms.count_stats.queue_stats import QueueStats


class CountPipeline:
    fivepseq_counts = None
    fivepseq_out = None
    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

    def __init__(self, fivepseq_counts, fivepseq_out):
        self.fivepseq_counts = fivepseq_counts
        self.fivepseq_out = fivepseq_out

    def run(self):

        # info
        self.logger.info("\n\nFivepseq started for bam file %s\n\n"
                         % os.path.basename(self.fivepseq_counts.alignment.bam_file))

        self.run_transcript_assambly()
        cancel = self.run_transcript_and_count_descriptors()

        if not cancel:
            self.run_outliers()
            self.run_start_stop_dicts()
            self.run_full_length_counts()
            self.run_count_stats()
            self.run_meta_counts()
            self.run_frame_counts()
            self.run_transcript_indices()
            self.run_codon_counts()
            self.run_aa_island_counts()
            self.run_codon_stats()
            self.run_loci_counts()
            self.run_queue_analysis()
            self.run_out_of_frame_codon_stats()

        # sanity check
        success = self.fivepseq_out.sanity_check_for_counts()
        if success:
            self.logger.info(
                "\n\nFivepseq successfully finished for bam file %s\n\n"
                % os.path.basename(self.fivepseq_counts.alignment.bam_file))
        else:
            self.logger.info(
                "\n\nFivepseq finished for bam file %s.\n Some files failed to be generated. Check those in %s.\n\n"
                % (os.path.basename(self.fivepseq_counts.alignment.bam_file),
                   self.fivepseq_out.get_file_path(FivePSeqOut.FAILED_COUNT_FILES_LIST)))

    def run_transcript_assambly(self):
        #   transcript assembly
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE)):
            self.fivepseq_out.write_transcript_assembly_to_file(
                self.fivepseq_counts.annotation.get_transcript_assembly(span_size=0),
                self.fivepseq_out.TRANSCRIPT_ASSEMBLY_FILE,
                self.fivepseq_counts,
                span_size=0
            )

    def run_transcript_and_count_descriptors(self):
        cancel = False
        #   gene sets: always re-write gene sets, because those can change at each run
        # if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.GENESET_FILE)):
        if self.fivepseq_counts.annotation.gs_transcript_dict is not None:
            self.fivepseq_out.write_geneset_transcript_dict_to_file(self.fivepseq_counts.annotation.gs_transcript_dict,
                                                                    self.fivepseq_out.GENESET_FILE)

        # case 1 transcript descriptor are not there: generate
        # case 2 transcript descriptors are there, but both count_distribution and outliers are not there: generate
        # case 3 transcript descriptors are there, count_distriubtion is there: set count_distribution
        # case 4 transcript descriptors are there, outliers are there: set ouli

        # transcript descriptors
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)) or \
                not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_DISTRIBUTION_FILE)):

            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_transcript_descriptors(),
                                               self.fivepseq_out.TRANSCRIPT_DESCRIPTORS_FILE)
            self.fivepseq_out.write_dict(self.fivepseq_counts.get_count_distribution_dict(),
                                         self.fivepseq_out.COUNT_DISTRIBUTION_FILE)

        else:
            count_distribution = CountManager.read_count_dict(self.fivepseq_out.get_file_path(
                self.fivepseq_out.COUNT_DISTRIBUTION_FILE))
            self.fivepseq_counts.set_count_distribution_dict(count_distribution)

        count_distribution = self.fivepseq_counts.get_count_distribution()
        if len(count_distribution) == 0:
            self.logger.warning(
                "No reads found in coding regions. Fivepseq will skip the rest of calculations.")
            cancel = True

        return cancel

    def run_outliers(self):
        self.fivepseq_out.write_dict(self.fivepseq_counts.get_count_distribution_dict(),
                                     self.fivepseq_out.COUNT_DISTRIBUTION_FILE)
        self.fivepseq_out.write_vector([self.fivepseq_counts.get_outlier_lower()],
                                       self.fivepseq_out.OUTLIER_LOWER_FILE)

    def run_start_stop_dicts(self):
        #   start codon dictionary
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.START_CODON_DICT_FILE)):
            self.fivepseq_out.write_dict(self.fivepseq_counts.get_start_codon_dict(),
                                         self.fivepseq_out.START_CODON_DICT_FILE)

        #   stop codon dictionary
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.TERM_CODON_DICT_FILE)):
            self.fivepseq_out.write_dict(self.fivepseq_counts.get_stop_codon_dict(),
                                         self.fivepseq_out.TERM_CODON_DICT_FILE)

    def run_full_length_counts(self):
        #   load or generate full length counts
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_FULL_FILE)):
            self.fivepseq_out.write_vector_list(
                self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.FULL_LENGTH),
                self.fivepseq_out.COUNT_FULL_FILE)

        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.OUTLIERS_DF)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_outliers_df(),
                                               self.fivepseq_out.OUTLIERS_DF)

    def run_count_stats(self):
        # count stats
        count_stats = CountStats(self.fivepseq_counts, self.fivepseq_out, config)
        count_stats.count_stats()

    def run_meta_counts(self):
        #   terminal counts
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_TERM_FILE)):
            term_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.TERM)
            self.fivepseq_out.write_vector_list(term_counts, self.fivepseq_out.COUNT_TERM_FILE)

        #   start counts
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.COUNT_START_FILE)):
            start_counts = self.fivepseq_counts.get_count_vector_list(FivePSeqCounts.START)
            self.fivepseq_out.write_vector_list(start_counts, self.fivepseq_out.COUNT_START_FILE)

        #   meta counts term
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.META_COUNT_TERM_FILE)):
            self.fivepseq_out.write_series_to_file(self.fivepseq_counts.get_meta_count_series(FivePSeqCounts.TERM),
                                                   self.fivepseq_out.META_COUNT_TERM_FILE)

        #   meta counts start
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.META_COUNT_START_FILE)):
            self.fivepseq_out.write_series_to_file(self.fivepseq_counts.get_meta_count_series(FivePSeqCounts.START),
                                                   self.fivepseq_out.META_COUNT_START_FILE)

    def run_frame_counts(self):
        #   frame counts start
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_START_FILE)):
            frame_count_df_start = self.fivepseq_counts.get_frame_counts_df(FivePSeqCounts.START)
            self.fivepseq_out.write_df_to_file(frame_count_df_start, self.fivepseq_out.FRAME_COUNTS_START_FILE)

        #   frame counts term
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.FRAME_COUNTS_TERM_FILE)):
            frame_count_df_term = self.fivepseq_counts.get_frame_counts_df(FivePSeqCounts.TERM)
            self.fivepseq_out.write_df_to_file(frame_count_df_term, self.fivepseq_out.FRAME_COUNTS_TERM_FILE)

    def run_transcript_indices(self):
        #   canonical transcript indices
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)):
            self.fivepseq_out.write_vector(self.fivepseq_counts.canonical_transcript_index,
                                           self.fivepseq_out.CANONICAL_TRANSCRIPT_INDEX_FILE)

    def run_codon_counts(self):
        #   amino acid pauses
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.AMINO_ACID_PAUSES_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_amino_acid_pauses(),
                                               # generate more than needed for visualization
                                               self.fivepseq_out.AMINO_ACID_PAUSES_FILE)

        #   codon pauses

        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.CODON_PAUSES_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_codon_pauses(),
                                               self.fivepseq_out.CODON_PAUSES_FILE)

        #   dicodon pauses

        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.DICODON_PAUSES_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_dicodon_pauses(),
                                               self.fivepseq_out.DICODON_PAUSES_FILE)
        #   dipeptide pauses

        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.DIPEPTIDE_PAUSES_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_dipeptide_pauses(),
                                               self.fivepseq_out.DIPEPTIDE_PAUSES_FILE)

        #   tricodon pauses

        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.TRICODON_PAUSES_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_tricodon_pauses(),
                                               self.fivepseq_out.TRICODON_PAUSES_FILE)
        #   tripeptide pauses

        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.TRIPEPTIDE_PAUSES_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_tripeptide_pauses(),
                                               self.fivepseq_out.TRIPEPTIDE_PAUSES_FILE)

    def run_aa_island_counts(self):

        if "aa_islands" in config.args and config.args.aa_islands is True:

            if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.AMINO_ACID_ISLAND_PAUSES_FILE)):
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_amino_acid_island_pauses(),
                                                   # generate more than needed for visualization
                                                   self.fivepseq_out.AMINO_ACID_ISLAND_PAUSES_FILE)


    def run_codon_stats(self):
        #   amino acid stats
        if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.AMINO_ACID_STATS_FILE)):
            self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_amino_acid_stats(),
                                               # generate more than needed for visualization
                                               self.fivepseq_out.AMINO_ACID_STATS_FILE)

    def run_out_of_frame_codon_stats(self):

        if "oof" in config.args and config.args.oof is True:
            #   frame +1
            self.logger.info("Computing out of frame codon counts (F2)")
            if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.CODON_PAUSES_F2_FILE)):
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_out_of_frame_codon_pauses (fstart=1),
                                                   self.fivepseq_out.CODON_PAUSES_F2_FILE)
            self.logger.info("Computing out of frame codon counts (F0)")
            # frame -1
            if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.CODON_PAUSES_F0_FILE)):
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_out_of_frame_codon_pauses(fstart=2),
                                                   self.fivepseq_out.CODON_PAUSES_F0_FILE)


    def run_loci_counts(self):
        #   loci pauses
        read_locations = FivePSeqCounts.READ_LOCATIONS_ALL
        loci_pauses_filename = self.fivepseq_out.LOCI_PAUSES_FILE_PREFIX + read_locations + ".txt"
        if not self.skip(self.fivepseq_out.get_file_path(loci_pauses_filename)):
            if self.fivepseq_counts.loci_file is not None:
                self.fivepseq_out.write_series_to_file(
                    self.fivepseq_counts.get_pauses_from_loci(self.fivepseq_counts.loci_file,
                                                              read_locations=read_locations),
                    loci_pauses_filename)
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_loci_overlaps_df(),
                                                   self.fivepseq_out.LOCI_OVERLAPS_FILE_PREFIX + read_locations)

        read_locations = FivePSeqCounts.READ_LOCATIONS_3UTR
        loci_pauses_filename = self.fivepseq_out.LOCI_PAUSES_FILE_PREFIX + read_locations + ".txt"
        if not self.skip(self.fivepseq_out.get_file_path(loci_pauses_filename)):
            if self.fivepseq_counts.loci_file is not None:
                self.fivepseq_out.write_series_to_file(
                    self.fivepseq_counts.get_pauses_from_loci(self.fivepseq_counts.loci_file,
                                                              read_locations=read_locations),
                    loci_pauses_filename
                )
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_loci_overlaps_df(),
                                                   self.fivepseq_out.LOCI_OVERLAPS_FILE_PREFIX + read_locations)

        read_locations = FivePSeqCounts.READ_LOCATIONS_5UTR
        loci_pauses_filename = self.fivepseq_out.LOCI_PAUSES_FILE_PREFIX + read_locations + ".txt"
        if not self.skip(self.fivepseq_out.get_file_path(loci_pauses_filename)):
            if self.fivepseq_counts.loci_file is not None:
                self.fivepseq_out.write_series_to_file(
                    self.fivepseq_counts.get_pauses_from_loci(self.fivepseq_counts.loci_file,
                                                              read_locations=read_locations),
                    loci_pauses_filename)
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_loci_overlaps_df(),
                                                   self.fivepseq_out.LOCI_OVERLAPS_FILE_PREFIX + read_locations)

        read_locations = FivePSeqCounts.READ_LOCATIONS_CDS
        loci_pauses_filename = self.fivepseq_out.LOCI_PAUSES_FILE_PREFIX + read_locations + ".txt"
        if not self.skip(self.fivepseq_out.get_file_path(loci_pauses_filename)):
            if self.fivepseq_counts.loci_file is not None:
                self.fivepseq_out.write_series_to_file(
                    self.fivepseq_counts.get_pauses_from_loci(self.fivepseq_counts.loci_file,
                                                              read_locations=read_locations),
                    loci_pauses_filename
                )
                self.fivepseq_out.write_df_to_file(self.fivepseq_counts.get_loci_overlaps_df(),
                                                   self.fivepseq_out.LOCI_OVERLAPS_FILE_PREFIX + read_locations)

    def run_queue_analysis(self):
        if "queue" in config.args and config.args.queue is not None:
            if not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.QUEUE_STATS_FILE))\
                    and not self.skip(self.fivepseq_out.get_file_path(self.fivepseq_out.OUT_OF_FRAME_STATS_FILE)):
                queue_stats = QueueStats(self.fivepseq_counts, self.fivepseq_out, config)
                try:
                    period = int(config.args.queue)
                except:
                    logging.warning("Illegal argument %s to queue option" % config.args.queue)
                    return
                rq_df = queue_stats.store_periods(period)
                self.fivepseq_out.write_df_to_file(rq_df, self.fivepseq_out.get_file_path(
                    self.fivepseq_out.QUEUE_STATS_FILE))
                out_of_frame_df = self.fivepseq_counts.store_out_of_frame_stats(period)
                self.fivepseq_out.write_df_to_file(out_of_frame_df, self.fivepseq_out.get_file_path(
                    self.fivepseq_out.OUT_OF_FRAME_STATS_FILE))

    def skip(self, file):
        if (self.fivepseq_out.conflict_mode == config.ADD_FILES) & (
                os.path.exists(self.fivepseq_out.get_file_path(file))):
            self.logger.info(
                "Skipping existing file %s" % file)
            return True

        return False
