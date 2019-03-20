import colorlover as cl
import glob
import logging
import os
import pandas as pd
from bokeh.io import export_svgs, export_png
from bokeh.plotting import figure

import config
from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from fivepseq.util.writers import FivePSeqOut
from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_heatmap_grid, bokeh_frame_barplots, \
    bokeh_composite, bokeh_fft_plot


class VizPipeline:
    FILTER_TOP_POPULATED = "populated"
    FILTER_CANONICAL_TRANSCRIPTS = "canonical"

    METACOUNTS_TERM = "_metacounts_term"
    METACOUNTS_START = "_metacounts_start"
    TRIANGLE_TERM = "_triangle_term"
    TRIANGLE_START = "_triangle_start"
    FRAME_TERM = "_frames_term"
    FRAME_START = "_frames_start"
    AMINO_ACID_PAUSES = "_amino_acid_pauses"
    AMINO_ACID_PAUSES_SCALED = "_amino_acid_pauses_scaled"
    FFT_TERM = "_fft_term"
    FFT_START = "_fft_start"




    png_dir = "png"
    svg_dir = "svg"

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

    phantomjs_installed = None

    scatter_term = None
    p_scatter_start = None
    p_triangle_term = None
    p_triangle_start = None
    p_aa_heatmap = None
    p_aa_heatmap_scaled = None
    p_frame_barplots_term = None
    p_frame_barplots_start = None
    p_loci_meta_counts = None
    p_fft_plot_start = None
    p_fft_plot_term = None

    p_scatter_term_combined = None
    p_scatter_start_combined = None
    p_triangle_term_combined = None
    p_triangle_start_combined = None
    p_aa_heatmaps_combined = None

    def __init__(self, args):
        self.args = args

    def run(self):
        try:
            self.prepare_output_folders()
        except Exception as e:
            err_msg = "Problem with plotting: could not create folders for exporting images: %s" % str(e)
            logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).error(err_msg)
            raise Exception(err_msg)

        if hasattr(self.args, 'sd') and  self.args.sd is not None:
            self.logger.info("Plotting single sample: %s" % self.args.sd)

            self.initialize_data()
            self.logger.info("Finished reading data counts.")

            if self.args.t is None:
                title = os.path.basename(self.args.sd)
            else:
                title = self.args.t

            self.plot_single_sample(title, self.args.o)

        else:
            self.logger.info("Plotting multiple samples:")
            for d in glob.glob(self.args.md):
                if os.path.isdir(d):
                    self.logger.info("\t%s" % d)

            self.initialize_data()
            self.logger.info("Finished reading data counts.")
            combined = True
            if self.args.tf is None:
                self.combine_counts()
            else:
                combined = False
            self.plot_multiple_samples(combined)

    def prepare_output_folders(self):
        if not os.path.exists(self.args.o):
            try:
                os.mkdir(self.args.o)
            except Exception as e:
                raise Exception("Output directory %s could not be created: %s" % (self.args.o, str(e)))

        if self.is_phantomjs_installed():
            self.png_dir = os.path.join(self.args.o, "png")
            if not os.path.exists(self.png_dir):
                try:
                    os.mkdir(self.png_dir)
                except Exception as e:
                    raise Exception("Output directory %s could not be created: %s" % (self.png_dir, str(e)))

            self.svg_dir = os.path.join(self.args.o, "svg")
            if not os.path.exists(self.svg_dir):
                try:
                    os.mkdir(self.svg_dir)
                except Exception as e:
                    raise Exception("Output directory %s could not be created: %s" % (self.svg_dir, str(e)))

    def initialize_data(self):
        if hasattr(self.args, 'sd') and self.args.sd is not None:
            if self.args.t is None:
                sample = os.path.basename(self.args.sd)
            else:
                sample = self.args.t
            self.samples.append(sample)
            self.update_dicts(sample, self.args.sd)
        else:
            for d in glob.glob(self.args.md):
                if os.path.isdir(d):
                    if d[-1] == "/":
                        d = d[0:len(d) - 1]
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

    def plot_single_sample(self, title, plot_dir):
        self.make_single_sample_plots(title)

        # TODO skip tables for now as those are not properly drawn
        bokeh_composite(title, [self.p_scatter_start, self.p_scatter_term,
                                self.p_triangle_term, self.p_triangle_start,
                                self.p_aa_heatmap, None,
                                self.p_aa_heatmap_scaled, None,
                                self.p_frame_barplots_term, self.p_frame_barplots_start,
                                self.p_loci_meta_counts, None,
                                self.p_fft_plot_start, self.p_fft_plot_term],
                        os.path.join(plot_dir, title + ".html"), 2)

    def plot_multiple_samples(self, combined=True):
        self.logger.info("Generating plots for single samples")
        if self.args.t is None:
            self.args.t = os.path.basename(os.path.dirname(self.args.md)) + "_" + os.path.basename(self.args.md)

        self.make_single_sample_plots(self.args.t)
        if combined:
            self.make_combined_plots(self.args.t)

        if combined:

            # TODO skip tables for now as those are not properly drawn
            bokeh_composite(self.args.t,
                            [self.p_scatter_start, self.p_scatter_term,
                             self.p_triangle_term, self.p_triangle_start,
                             self.p_aa_heatmap, None, None, None,
                             self.p_aa_heatmap_scaled, None, None, None,
                             self.p_frame_barplots_term, None, None, None,
                             self.p_frame_barplots_start, None, None, None,
                             self.p_fft_plot_start, self.p_fft_plot_term, None, None,
                             self.p_scatter_start_combined, self.p_scatter_term_combined,
                             self.p_triangle_term_combined, self.p_triangle_start_combined,
                             self.p_aa_heatmaps_combined, None, None, None],
                            os.path.join(self.args.o, self.args.t + ".html"), 4)
        else:
            bokeh_composite(self.args.t,
                            [self.p_scatter_start, self.p_scatter_term,
                             self.p_triangle_term, self.p_triangle_start,
                             self.p_aa_heatmap, None, None, None,
                             self.p_aa_heatmap_scaled, None, None, None,
                             self.p_frame_barplots_term, None, None, None,
                             self.p_frame_barplots_start, None, None, None,
                             self.p_fft_plot_start, self.p_fft_plot_term, None, None],
                            os.path.join(self.args.o, self.args.t + ".html"), 4)

    def is_phantomjs_installed(self):
        if self.phantomjs_installed is not None:
            return self.phantomjs_installed

        # prepare some data
        x = [1, 2, 3, 4, 5]
        y = [6, 7, 2, 4, 5]

        # output to static HTML file
        # output_file("lines.html")

        # create a new plot with a title and axis labels
        p = figure(title="simple line example", x_axis_label='x', y_axis_label='y')

        # add a line renderer with legend and line thickness
        p.line(x, y, legend="Temp.", line_width=2)

        p.output_backend = "svg"
        try:
            test_file_name = "fivepseq.phantom.test.svg"
            export_svgs(p, filename=test_file_name)
            self.phantomjs_installed = True
            os.remove(test_file_name)
        except:
            self.phantomjs_installed = False
            # TODO in the future fivepseq should attempt to install phantomjs from the very beginning
            logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).warning(
                "It seems like phantomjs is not installed no your system. "
                "Files may not be exported in svg and png formats, while html will still be available for viewing."
                "To install phantomjs, run 'conda install phantomjs selenium pillow'")
        return self.phantomjs_installed

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

    def make_single_sample_plots(self, title):
        self.p_scatter_term = bokeh_scatter_plot(title + self.METACOUNTS_TERM, FivePSeqCounts.TERM,
                                                 self.meta_count_term_dict,
                                                 self.colors_dict,
                                                 png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_scatter_start = bokeh_scatter_plot(title + self.METACOUNTS_START, FivePSeqCounts.START,
                                                  self.meta_count_start_dict,
                                                  self.colors_dict,
                                                  png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_triangle_term = bokeh_triangle_plot(title + self.TRIANGLE_TERM,
                                                   self.frame_count_term_dict,
                                                   self.colors_dict,
                                                   png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_triangle_start = bokeh_triangle_plot(title + self.TRIANGLE_START,
                                                    self.frame_count_start_dict,
                                                    self.colors_dict,
                                                    png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_aa_heatmap = bokeh_heatmap_grid(title + self.AMINO_ACID_PAUSES,
                                               self.amino_acid_df_dict,
                                               png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_aa_heatmap_scaled = bokeh_heatmap_grid(title + self.AMINO_ACID_PAUSES_SCALED,
                                                      self.amino_acid_df_dict, scale=True,
                                                      png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_frame_barplots_term = bokeh_frame_barplots(title + self.FRAME_TERM,
                                                          self.frame_count_term_dict,
                                                          self.frame_stats_df_dict,
                                                          self.colors_dict,
                                                          png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_frame_barplots_start = bokeh_frame_barplots(title + self.FRAME_START,
                                                           self.frame_count_start_dict,
                                                           self.frame_stats_df_dict,
                                                           self.colors_dict,
                                                           png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_loci_meta_counts = bokeh_scatter_plot(title + "_loci_meta_counts",
                                                     "loci",
                                                     self.loci_meta_counts_dict,
                                                     self.colors_dict,
                                                     png_dir=self.png_dir, svg_dir=self.svg_dir)

        # fft plots
        self.p_fft_plot_start = bokeh_fft_plot(title + self.FFT_START, "start",
                                               self.fft_signal_start_dict,
                                               self.colors_dict,
                                               png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_fft_plot_term = bokeh_fft_plot(title + self.FFT_TERM, "term",
                                              self.fft_signal_term_dict,
                                              self.colors_dict,
                                              png_dir=self.png_dir, svg_dir=self.svg_dir)

    def make_combined_plots(self, title):
        # Combined plots
        self.p_scatter_term_combined = bokeh_scatter_plot(title + self.METACOUNTS_TERM + "_combined", FivePSeqCounts.TERM,
                                                          {"combined": self.meta_count_term_combined},
                                                          self.combined_color_dict,
                                                          png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_scatter_start_combined = bokeh_scatter_plot(title + self.METACOUNTS_START + "_combined", FivePSeqCounts.START,
                                                           {self.COMBINED: self.meta_count_start_combined},
                                                           self.combined_color_dict,
                                                           png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_triangle_term_combined = bokeh_triangle_plot(title + self.TRIANGLE_TERM + "_combined",
                                                            {self.COMBINED: self.frame_count_TERM_combined},
                                                            self.combined_color_dict,
                                                            png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_triangle_start_combined = bokeh_triangle_plot(title + self.TRIANGLE_START + "_combined",
                                                             {self.COMBINED: self.frame_count_START_combined},
                                                             self.combined_color_dict,
                                                             png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_aa_heatmaps_combined = bokeh_heatmap_grid(title + self.AMINO_ACID_PAUSES + "_combined",
                                                         {self.COMBINED: self.amino_acid_df_combined},
                                                         png_dir=self.png_dir, svg_dir=self.svg_dir)
