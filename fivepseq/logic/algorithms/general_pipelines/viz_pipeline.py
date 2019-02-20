import glob
import os

import pandas as pd
import colorlover as cl
from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_heatmap_grid, bokeh_frame_barplots, \
    bokeh_composite

import config
from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from fivepseq.util.writers import FivePSeqOut
from fivepseq.viz.visualization_pipeline import plot_single_sample


class VizPipeline:
    fivepseq_counts = None
    fivepseq_out = None
    args = None

    samples = []

    transcript_count_full_dict = {}
    meta_count_term_dict = {}
    meta_count_start_dict = {}
    amino_acid_df_dict = {}
    frame_count_dict = {}
    frame_count_dict_START = {}

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

    def run(self, logger, filter=None):
        if self.args.sd is not None:
            logger.info("Plotting single sample: %s" % self.args.sd)
            if self.args.t is None:
                title = os.path.basename(self.args.sd)
            else:
                title = self.args.t
            plot_single_sample(title, title + "_plot",
                               self.args.sd, self.args.o, self.args.span,
                               filter)

        else:
            logger.info("Plotting multiple samples:")
            for d in glob.glob(self.args.md):
                logger.info("\t%s" % d)

            self.initialize_dicts(logger)
            logger.info("Finished reading data counts.")

            self.combine_counts(logger)
            self.plot_multiple_samples(logger)

    def initialize_dicts(self, logger):
        for d in glob.glob(self.args.md):
            sample = os.path.basename(d)
            self.samples.append(sample)
            logger.info("reading counts for sample: %s" % sample)

            file = os.path.join(d, FivePSeqOut.COUNT_FULL_FILE)
            self.transcript_count_full_dict.update(
                {sample: CountManager.read_counts_as_list(file)})

            file = os.path.join(d, FivePSeqOut.META_COUNT_TERM_FILE)
            self.meta_count_term_dict.update(
                {sample: pd.read_csv(file, sep="\t", header=None, names=["D", "C"])})

            file = os.path.join(d, FivePSeqOut.META_COUNT_START_FILE)
            self.meta_count_start_dict.update(
                {sample: pd.read_csv(file, sep="\t", header=None, names=["D", "C"])})

            file = os.path.join(d, FivePSeqOut.FRAME_COUNTS_TERM_FILE)
            self.frame_count_dict.update(
                {sample: pd.read_csv(file, sep="\t")})

            file = os.path.join(d, FivePSeqOut.FRAME_COUNTS_START_FILE)
            self.frame_count_dict_START.update(
                {sample: pd.read_csv(file, sep="\t")})

            file = os.path.join(d, FivePSeqOut.AMINO_ACID_PAUSES_FILE)
            if os.path.exists(file):
                self.amino_acid_df_dict.update(
                    {sample: pd.read_csv(file, sep="\t", header=0, index_col=0)})
            else:
                logger.warning("Could not find amino acid pauses file %s. No plots will be generated for it." % file)
                self.amino_acid_df_dict.update({sample: None})

            self.colors_dict = dict(
                zip(self.transcript_count_full_dict.keys(),
                    self.large_colors_list[0:len(self.samples)]))

    def combine_counts(self, logger):
        logger.info("Combining data counts.")
        for key in self.meta_count_start_dict.keys():
            df_start = self.meta_count_start_dict.get(key)
            df_term = self.meta_count_term_dict.get(key)
            frame_start = self.frame_count_dict_START.get(key)
            frame_term = self.frame_count_dict.get(key)
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

    def plot_multiple_samples(self, logger):
        logger.info("Generating plots for single samples")
        p_scatter_term = bokeh_scatter_plot("meta-counts_TERM", FivePSeqCounts.TERM,
                                            self.meta_count_term_dict, self.colors_dict)

        p_scatter_start = bokeh_scatter_plot("meta-counts_START", FivePSeqCounts.START,
                                             self.meta_count_start_dict, self.colors_dict)

        p_triangle = bokeh_triangle_plot("frame_preferences_TERM",
                                         self.frame_count_dict, self.colors_dict)

        p_triangle_START = bokeh_triangle_plot("frame_preferences_START",
                                               self.frame_count_dict_START, self.colors_dict)

        p_aa_heatmaps = bokeh_heatmap_grid("amino_acid_pauses", self.amino_acid_df_dict)

        p_frame_barplots = bokeh_frame_barplots("frame_histograms_TERM",
                                                self.frame_count_dict, self.colors_dict)

        p_frame_barplots_start = bokeh_frame_barplots("frame_histograms_START",
                                                      self.frame_count_dict_START,
                                                      self.colors_dict)

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

        p_frame_barplots_combined = bokeh_frame_barplots("frame_histograms_TERM_combined",
                                                         {self.COMBINED: self.frame_count_TERM_combined},
                                                         self.combined_color_dict)
        p_frame_barplots_start_combined = bokeh_frame_barplots("frame_histograms_START",
                                                               {self.COMBINED: self.frame_count_START_combined},
                                                               self.combined_color_dict)

        if self.args.t is None:
            self.args.t = os.path.basename(os.path.dirname(self.args.md)) + "_" + os.path.basename(self.args.md)

        bokeh_composite(self.args.t,
                        [p_scatter_start, p_scatter_term,
                         p_triangle, p_triangle_START,
                         p_aa_heatmaps, None, None, None,
                         p_frame_barplots, None, None, None,
                         p_frame_barplots_start, None, None, None,
                         p_scatter_start_combined, p_scatter_term_combined,
                         p_triangle_combined, p_triangle_START_combined,
                         p_aa_heatmaps_combined, None, None, None,
                         p_frame_barplots_combined, None, None, None,
                         p_frame_barplots_start_combined, None, None, None],
                        os.path.join(self.args.o, self.args.t + ".html"), 4)
