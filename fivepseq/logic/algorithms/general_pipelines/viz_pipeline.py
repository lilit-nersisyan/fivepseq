import glob
import os

import logging
import pandas as pd
import colorlover as cl
from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_heatmap_grid, bokeh_frame_barplots, \
    bokeh_composite, bokeh_fft_plot, bokeh_table

import config
from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from fivepseq.util.writers import FivePSeqOut


class VizPipeline:
    FILTER_TOP_POPULATED = "populated"
    FILTER_CANONICAL_TRANSCRIPTS = "canonical"

    logger = logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER)

    fivepseq_counts = None
    args = None

    samples = []

    meta_count_term_dict = {}
    meta_count_start_dict = {}
    count_vector_list_start_dict = {}
    count_vector_list_term_dict = {}
    amino_acid_df_dict = {}
    frame_count_term_dict = {}
    frame_count_start_dict = {}
    frame_stats_df_dict = {}
    fft_signal_start_dict = {}
    fft_signal_term_dict = {}
    loci_meta_counts_dict = {}
    data_summary_dict = {}
    transcript_index = None

    COMBINED = "combined"
    meta_count_start_combined = pd.DataFrame()
    meta_count_term_combined = pd.DataFrame()
    frame_count_START_combined = pd.DataFrame()
    frame_count_TERM_combined = pd.DataFrame()
    amino_acid_df_combined = pd.DataFrame()

    large_colors_list = (cl.to_numeric(cl.scales['8']['qual']['Paired']) +
                         cl.to_numeric(cl.scales['8']['qual']['Set1']) +
                         cl.to_numeric(cl.scales['5']['qual']['Set3']))
    colors_dict = None
    combined_color_dict = {COMBINED: cl.to_numeric(cl.scales['9']['qual']['Set3'])[3]}

    def __init__(self, args):
        self.args = args

    def run(self):
        if self.args.sd is not None:
            self.logger.info("Plotting single sample: %s" % self.args.sd)

            self.initialize_data()
            self.logger.info("Finished reading data counts.")

            if self.args.t is None:
                title = os.path.basename(self.args.sd)
            else:
                title = self.args.t

            self.plot_single_sample(title, title, self.args.o)

        else:
            self.logger.info("Plotting multiple samples:")
            for d in glob.glob(self.args.md):
                self.logger.info("\t%s" % d)

            self.initialize_data()
            self.logger.info("Finished reading data counts.")
            combined = True
            if self.args.tf is None:
                self.combine_counts()
            else:
                combined = False
            self.plot_multiple_samples(combined)

    def initialize_data(self):
        if self.args.sd is not None:
            if self.args.t is None:
                sample = os.path.basename(self.args.sd)
            else:
                sample = self.args.t
            self.samples.append(sample)
            self.update_dicts(sample, self.args.sd)
        else:
            for d in glob.glob(self.args.md):
                if d[-1] == "/":
                    d = d[0:len(d)-1]
                sample = os.path.basename(d)
                logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).debug("Sample: %s" % sample)
                self.samples.append(sample)
                self.update_dicts(sample, d)

        self.colors_dict = dict(
            zip(self.meta_count_start_dict.keys(),
                self.large_colors_list[0:len(self.samples)]))

    def update_dicts(self, sample, directory):
        self.logger.info("reading counts for sample: %s" % sample)
        fivepseq_out = FivePSeqOut(directory)

        self.data_summary_dict.update({sample: self.read_data_summary(fivepseq_out)})
        self.meta_count_start_dict.update({sample: self.read_meta_count_start(fivepseq_out)})
        self.meta_count_term_dict.update({sample: self.read_meta_count_term(fivepseq_out)})

        self.frame_count_term_dict.update({sample: self.read_frame_count_term(fivepseq_out)})
        self.frame_count_start_dict.update({sample: self.read_frame_count_start(fivepseq_out)})
        self.frame_stats_df_dict.update({sample: self.read_frame_stats_df(fivepseq_out)})
        self.amino_acid_df_dict.update({sample: self.read_amino_acid_df(fivepseq_out)})

        self.fft_signal_start_dict.update({sample: self.read_fft_signal_start(fivepseq_out)})
        self.fft_signal_term_dict.update({sample: self.read_fft_signal_term(fivepseq_out)})

        self.count_vector_list_start_dict.update({sample: self.read_count_vector_list_start(fivepseq_out)})
        self.count_vector_list_term_dict.update({sample: self.read_count_vector_list_term(fivepseq_out)})

        self.loci_meta_counts_dict.update({sample: self.read_loci_meta_counts(fivepseq_out)})

        if self.args.tf is not None:
            filter = self.args.tf
            if filter == self.FILTER_TOP_POPULATED:
                self.logger.info("Applying filter %s" % filter)
                self.transcript_index = CountManager.top_populated_count_vector_indices(
                    self.count_vector_list_term_dict.get(sample), self.args.span, 10000)
            elif filter == self.FILTER_CANONICAL_TRANSCRIPTS:
                self.logger.info("Applying filter %s" % filter)
                self.transcript_index = CountManager.canonical_transcript_indices(directory)

        if self.transcript_index is not None:
            self.logger.info("Number of filtered transcripts: %d" % len(self.transcript_index))
            self.frame_count_term_dict[sample] = self.frame_count_term_dict[sample].iloc[self.transcript_index,]
            self.frame_count_start_dict[sample] = self.frame_count_start_dict[sample].iloc[self.transcript_index,]
            self.frame_stats_df_dict[sample] = None

            self.count_vector_list_term_dict[sample] = [self.count_vector_list_term_dict[sample][i]
                                                        for i in self.transcript_index]
            self.count_vector_list_start_dict[sample] = [self.count_vector_list_start_dict[sample][i]
                                                         for i in self.transcript_index]
            self.meta_count_term_dict[sample] = CountManager.count_vector_to_df(
                CountManager.compute_meta_counts(self.count_vector_list_term_dict[sample]),
                FivePSeqCounts.TERM, self.args.span)
            self.meta_count_start_dict[sample] = CountManager.count_vector_to_df(
                CountManager.compute_meta_counts(self.count_vector_list_start_dict[sample]),
                FivePSeqCounts.TERM, self.args.span)

            # TODO amino acids pauses not subsettable

    def combine_counts(self):
        self.logger.info("Combining data counts.")
        for key in self.meta_count_start_dict.keys():
            df_start = self.meta_count_start_dict.get(key)
            df_term = self.meta_count_term_dict.get(key)
            frame_start = self.frame_count_start_dict.get(key)
            frame_term = self.frame_count_term_dict.get(key)
            amino_acid_df = self.amino_acid_df_dict.get(key)
            if len(self.meta_count_start_combined) == 0:
                self.meta_count_start_combined = df_start.copy()
                self.meta_count_term_combined = df_term.copy()
                self.frame_count_START_combined = frame_start.copy()
                self.frame_count_TERM_combined = frame_term.copy()
                self.amino_acid_df_combined = amino_acid_df.copy()
            else:
                self.meta_count_start_combined.C += df_start.C
                self.meta_count_term_combined.C += df_term.C
                self.frame_count_START_combined.loc[:, ('F0', 'F1', 'F2')] += frame_start.loc[:, ('F0', 'F1', 'F2')]
                self.frame_count_TERM_combined.loc[:, ('F0', 'F1', 'F2')] += frame_term.loc[:, ('F0', 'F1', 'F2')]
                if amino_acid_df is not None:
                    self.amino_acid_df_combined += amino_acid_df

    def plot_single_sample(self, sample_name, title, plot_dir):

        #    p_data_summary_table = bokeh_table(sample_name + "_data_summary",
        #                                       self.data_summary_dict)
        #    p_frame_stats_table = bokeh_table(sample_name + "frame_stats",
        #                                      self.frame_stats_df_dict)

        p_scatter_term = bokeh_scatter_plot(sample_name + "_term", FivePSeqCounts.TERM,
                                            self.meta_count_term_dict,
                                            self.colors_dict)
        p_scatter_start = bokeh_scatter_plot(sample_name + "_start", FivePSeqCounts.START,
                                             self.meta_count_start_dict,
                                             self.colors_dict)
        p_triangle = bokeh_triangle_plot(sample_name + "_triangle",
                                         self.frame_count_term_dict,
                                         self.colors_dict)

        p_triangle_START = bokeh_triangle_plot(sample_name + "_triangle_START",
                                               self.frame_count_start_dict,
                                               self.colors_dict)

        p_aa_heatmap = bokeh_heatmap_grid(sample_name + "_amino_acid_pauses",
                                          self.amino_acid_df_dict)

        p_aa_heatmap_scaled = bokeh_heatmap_grid(sample_name + "_amino_acid_pauses_scaled",
                                                 self.amino_acid_df_dict, scale=True)

        p_frame_barplots_term = bokeh_frame_barplots(sample_name + "_frame_histograms",
                                                     self.frame_count_term_dict,
                                                     self.frame_stats_df_dict,
                                                     self.colors_dict)
        p_frame_barplots_start = bokeh_frame_barplots(sample_name + "_frame_START_histograms",
                                                      self.frame_count_start_dict,
                                                      self.frame_stats_df_dict,
                                                      self.colors_dict)

        p_loci_meta_counts = bokeh_scatter_plot(sample_name + "_loci_meta_counts",
                                                "loci",
                                                self.loci_meta_counts_dict,
                                                self.colors_dict)

        # fft plots
        p_fft_plot_start = bokeh_fft_plot(sample_name + "_fft_start", "start",
                                          self.fft_signal_start_dict,
                                          self.colors_dict)

        p_fft_plot_term = bokeh_fft_plot(sample_name + "_fft_term", "term",
                                         self.fft_signal_term_dict,
                                         self.colors_dict)

        # TODO skip tables for now as those are not properly drawn
        bokeh_composite(title, [p_scatter_start, p_scatter_term,
                                p_triangle, p_triangle_START,
                                p_aa_heatmap, None,
                                p_aa_heatmap_scaled, None,
                                p_frame_barplots_term, p_frame_barplots_start,
                                p_loci_meta_counts, None,
                                p_fft_plot_start, p_fft_plot_term],
                        os.path.join(plot_dir, title + ".html"), 2)

    def plot_multiple_samples(self, combined = True):
        self.logger.info("Generating plots for single samples")
        if self.args.t is None:
            self.args.t = os.path.basename(os.path.dirname(self.args.md)) + "_" + os.path.basename(self.args.md)

        #    p_data_summary_table = bokeh_table("data_summary",
        #                                       self.data_summary_dict)
        #    p_frame_stats_table = bokeh_table("frame_stats",
        #                                      self.frame_stats_df_dict)
        p_scatter_term = bokeh_scatter_plot("meta-counts_TERM", FivePSeqCounts.TERM,
                                            self.meta_count_term_dict, self.colors_dict)

        p_scatter_start = bokeh_scatter_plot("meta-counts_START", FivePSeqCounts.START,
                                             self.meta_count_start_dict, self.colors_dict)

        p_triangle = bokeh_triangle_plot("frame_preferences_TERM",
                                         self.frame_count_term_dict, self.colors_dict)

        p_triangle_START = bokeh_triangle_plot("frame_preferences_START",
                                               self.frame_count_start_dict, self.colors_dict)

        p_aa_heatmaps = bokeh_heatmap_grid("amino_acid_pauses", self.amino_acid_df_dict)

        p_aa_heatmap_scaled = bokeh_heatmap_grid("amino_acid_pauses_scaled",
                                                 self.amino_acid_df_dict, scale=True)

        p_frame_barplots = bokeh_frame_barplots("frame_histograms_TERM",
                                                self.frame_count_term_dict,
                                                self.frame_stats_df_dict,
                                                self.colors_dict)

        p_frame_barplots_start = bokeh_frame_barplots("frame_histograms_START",
                                                      self.frame_count_start_dict,
                                                      self.frame_stats_df_dict,
                                                      self.colors_dict)

        p_fft_plot_start = bokeh_fft_plot("FFT_start", "start",
                                          self.fft_signal_start_dict,
                                          self.colors_dict)

        p_fft_plot_term = bokeh_fft_plot("FFT_term", "term",
                                         self.fft_signal_term_dict,
                                         self.colors_dict)

        if combined:
            # Combined plots
            p_scatter_term_combined = bokeh_scatter_plot("meta-counts_TERM_combined", FivePSeqCounts.TERM,
                                                         {"combined": self.meta_count_term_combined},
                                                         self.combined_color_dict)

            p_scatter_start_combined = bokeh_scatter_plot("meta-counts_START_combined", FivePSeqCounts.START,
                                                          {self.COMBINED: self.meta_count_start_combined},
                                                          self.combined_color_dict)

            p_triangle_combined = bokeh_triangle_plot("frame_preferences_TERM_combined",
                                                      {self.COMBINED: self.frame_count_TERM_combined},
                                                      self.combined_color_dict)

            p_triangle_START_combined = bokeh_triangle_plot("frame_preferences_START_combined",
                                                            {self.COMBINED: self.frame_count_START_combined},
                                                            self.combined_color_dict)

            p_aa_heatmaps_combined = bokeh_heatmap_grid("amino_acid_pauses_combined",
                                                        {self.COMBINED: self.amino_acid_df_combined})

            # p_frame_barplots_combined = bokeh_frame_barplots("frame_histograms_TERM_combined",
            #                                                 {self.COMBINED: self.frame_count_TERM_combined},
            #                                                 {self.COMBINED: None},
            #                                                 self.combined_color_dict)
            # p_frame_barplots_start_combined = bokeh_frame_barplots("frame_histograms_START",
            #                                                       {self.COMBINED: self.frame_count_START_combined},
            #                                                       {self.COMBINED: None},
            #                                                       self.combined_color_dict)


            # TODO skip tables for now as those are not properly drawn
            bokeh_composite(self.args.t,
                            [p_scatter_start, p_scatter_term,
                             p_triangle, p_triangle_START,
                             p_aa_heatmaps, None, None, None,
                             p_aa_heatmap_scaled, None, None, None,
                             p_frame_barplots, None, None, None,
                             p_frame_barplots_start, None, None, None,
                             p_fft_plot_start, p_fft_plot_term, None, None,
                             p_scatter_start_combined, p_scatter_term_combined,
                             p_triangle_combined, p_triangle_START_combined,
                             p_aa_heatmaps_combined, None, None, None],
                            os.path.join(self.args.o, self.args.t + ".html"), 4)
        else:
            bokeh_composite(self.args.t,
                            [p_scatter_start, p_scatter_term,
                             p_triangle, p_triangle_START,
                             p_aa_heatmaps, None, None, None,
                             p_aa_heatmap_scaled, None, None, None,
                             p_frame_barplots, None, None, None,
                             p_frame_barplots_start, None, None, None,
                             p_fft_plot_start, p_fft_plot_term, None, None],
                            os.path.join(self.args.o, self.args.t + ".html"), 4)

    def read_meta_count_term(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.META_COUNT_TERM_FILE)
        try:
            meta_count_term = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            meta_count_term = None
        return meta_count_term

    def read_meta_count_start(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.META_COUNT_START_FILE)
        try:
            meta_count_start = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            meta_count_start = None
        return meta_count_start

    def read_amino_acid_df(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.AMINO_ACID_PAUSES_FILE)
        try:
            amino_acid_df = CountManager.read_amino_acid_df(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            amino_acid_df = None
        return amino_acid_df

    def read_frame_count_term(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.FRAME_COUNTS_TERM_FILE)
        try:
            frame_count_term = CountManager.read_frame_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            frame_count_term = None
        return frame_count_term

    def read_frame_count_start(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.FRAME_COUNTS_START_FILE)
        try:
            frame_count_start = CountManager.read_frame_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            frame_count_start = None
        return frame_count_start

    def read_frame_stats_df(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.FRAME_STATS_DF_FILE)
        if os.path.exists(file):
            frame_stats_df = pd.read_csv(file, sep="\t", header=0, index_col=0)
        else:
            frame_stats_df = None
        return frame_stats_df

    def read_count_vector_list_term(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.COUNT_TERM_FILE)
        try:
            count_vector_list_term = CountManager.read_counts_as_list(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            count_vector_list_term = None
        return count_vector_list_term

    def read_count_vector_list_start(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.COUNT_START_FILE)
        try:
            count_vector_list_start = CountManager.read_counts_as_list(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            count_vector_list_start = None
        return count_vector_list_start

    def read_loci_meta_counts(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.LOCI_PAUSES_FILE)
        loci_meta_counts = None
        if os.path.exists(file):
            self.logger.info("Loci count file found")
            loci_meta_counts = CountManager.read_meta_counts(file)
        return loci_meta_counts

    def read_fft_signal_start(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.FFT_SIGNALS_START)
        try:
            fft_signals_start = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            fft_signals_start = None
        return fft_signals_start

    def read_fft_signal_term(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.FFT_SIGNALS_TERM)
        try:
            fft_signals_term = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            fft_signals_term = None
        return fft_signals_term

    def read_data_summary(self, fivepseq_out):
        file = fivepseq_out.get_file_path(FivePSeqOut.DATA_SUMMARY_FILE)
        try:
            data_summary = pd.read_csv(file, sep="\t", header=None, index_col=0)
        except:
            self.logger.warn("The file %s was not found, table will not be generated" % str(file))
            data_summary = None
        return data_summary
