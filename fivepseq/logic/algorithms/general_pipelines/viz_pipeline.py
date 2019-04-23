import colorlover as cl
import glob
import logging
import os
import pandas as pd
from bokeh.io import export_svgs, export_png
from bokeh.plotting import figure

from fivepseq import config
from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from fivepseq.util.writers import FivePSeqOut
from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_heatmap_grid, bokeh_frame_barplots, \
    bokeh_composite, bokeh_fft_plot
import numpy as np

from logic.structures.codons import Codons


class VizPipeline:
    FILTER_TOP_POPULATED = "populated"
    FILTER_CANONICAL_TRANSCRIPTS = "canonical"

    METACOUNTS_TERM = "_metacounts_term"
    METACOUNTS_START = "_metacounts_start"
    METACOUNTS_TERM_SCALED = "_metacounts_term_scaled"
    METACOUNTS_START_SCALED = "_metacounts_start_scaled"
    TRIANGLE_TERM = "_triangle_term"
    TRIANGLE_START = "_triangle_start"
    FRAME_TERM = "_frames_term"
    FRAME_START = "_frames_start"
    AMINO_ACID_PAUSES = "_amino_acid_pauses"
    AMINO_ACID_PAUSES_SCALED = "_amino_acid_pauses_scaled"
    FFT_TERM = "_fft_term"
    FFT_START = "_fft_start"

    count_folders = None
    title = "fivepseq_plot_canvas"
    png_dir = None
    svg_dir = None
    supplement_dir = "supplement"
    supplement_png_dir = None
    supplement_svg_dir = None

    logger = logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER)

    fivepseq_counts = None
    args = None

    samples = []

    meta_count_term_dict = {}
    meta_count_start_dict = {}
    count_vector_list_start_dict = {}
    count_vector_list_term_dict = {}
    amino_acid_df_dict = {}
    amino_acid_df_full_dict = {}
    codon_df_dict = {}
    codon_basesorted_df_dict = {}
    frame_count_term_dict = {}
    frame_count_start_dict = {}
    frame_stats_df_dict = {}
    fft_signal_start_dict = {}
    fft_signal_term_dict = {}
    loci_meta_counts_dict = {}
    data_summary_dict = {}
    transcript_index = None

    combine = False  # do not combine plots until data counts are successfully combined

    COMBINED = "combined"
    meta_count_start_combined = pd.DataFrame()
    meta_count_term_combined = pd.DataFrame()
    frame_count_START_combined = pd.DataFrame()
    frame_count_TERM_combined = pd.DataFrame()
    amino_acid_df_combined = pd.DataFrame()
    amino_acid_df_full_combined = pd.DataFrame()
    codon_df_combined = pd.DataFrame()
    codon_basesorted_df_combined = pd.DataFrame()

    large_colors_list = (cl.to_numeric(cl.scales['8']['qual']['Paired'][0:6]) +
                         cl.to_numeric(cl.scales['8']['qual']['Set1'][3:5]))
    #                     cl.to_numeric(cl.scales['5']['qual']['Set3']))
    # large_colors_list = ("#771155", "#AA4488", "#114477", "#4477AA", "#117777", "#44AAAA",
    #                     "#777711", "#AAAA44", "#774411", "#AA7744", "#771122", "#AA4455")
    colors_dict = None
    combined_color_dict = {COMBINED: cl.to_numeric(cl.scales['9']['qual']['Set3'])[3]}

    phantomjs_installed = None

    p_scatter_start = None
    p_triangle_term = None
    p_scatter_start_scaled = None
    p_triangle_term_scaled = None
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
    p_scatter_term_scaled_combined = None
    p_scatter_start_scaled_combined = None
    p_triangle_term_combined = None
    p_triangle_start_combined = None
    p_aa_heatmaps_combined = None

    # the distance to be used for plotting amino-acid heatmaps
    dist_for_amino_acid_heatmaps = 20

    def __init__(self, args, count_folders=None):
        """
        Initialize vizualization pipeline with arguments contained in args.
        If count_folders are provided explicitely, those will be used instead of sd and md arguments.

        :param args:
        :param count_folders:
        """
        self.args = args
        self.count_folders = count_folders

    def run(self):
        try:
            self.prepare_output_folders()
        except Exception as e:
            err_msg = "Problem with plotting: could not create folders for exporting images: %s" % str(e)
            self.logger.error(err_msg)
            raise e

        try:
            self.process_count_folders()
        except Exception as e:
            err_msg = "Exception while processing input count folders: %s" % str(e)
            self.logger.error(err_msg)
            raise e

        try:
            self.setup_title()
        except Exception as e:
            err_msg = "Exception while setting up title for plot canvas: %s" % str(e)
            self.logger.error(err_msg)
            raise e

        try:
            self.initialize_data()
        except Exception as e:
            err_msg = "Exception while reading data: %s" % str(e)
            self.logger.error(err_msg)
            raise e

        try:
            self.plot_multiple_samples()
        except Exception as e:
            err_msg = "Exception while plotting data: %s" % str(e)
            self.logger.error(err_msg)
            raise e

        try:
            self.write_supplement()
        except Exception as e:
            err_msg = "Exception while writing supplements: %s" % str(e)
            self.logger.error(err_msg)
            raise e

    def process_count_folders(self):
        """
        Read input as sd or md (if count_folders is not provided from within fivepseq)
        and append the count_folders list in self.

        Set up output file title.

        :return:
        """
        if self.count_folders is None:
            if hasattr(self.args, 'sd') and self.args.sd is not None:
                if not os.path.isdir(self.args.sd):
                    err_msg = "Provided sd %s is not a directory" % self.args.sd
                    self.logger.error(err_msg)
                    raise Exception(err_msg)
                self.count_folders.append(self.args.sd)

            elif hasattr(self.args, 'md') and self.args.sd is not None:
                for d in glob.glob(self.args.md):
                    if os.path.isdir(d):
                        if d[-1] == "/":
                            d = d[0:len(d) - 1]
                        self.count_folders.append(d)

        if len(self.count_folders) == 0:
            err_msg = "No count folders provided as input"
            self.logger.error(err_msg)
            raise Exception(err_msg)

        if len(self.count_folders) > 8:
            self.logger.info("The number of samples exceeds 8 (found %d). Only the first 8 will be plotted" % len(
                self.count_folders))
            self.count_folders = self.count_folders[0:8]

        self.logger.info("The following folders will be used for plotting")
        for d in self.count_folders:
            if os.path.isdir(d):
                self.logger.info("\t%s" % d)

    def setup_title(self):
        if not hasattr(config.args, 't') or config.args.t is None:
            config.args.t = os.path.basename(os.path.dirname(config.args.b)) + "_" + os.path.basename(config.args.b)
            config.args.t.replace("_*", "")
            config.args.t.replace("*", "")

        self.title = config.args.t

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

        self.supplement_dir = os.path.join(self.args.o, "supplement")
        if not os.path.exists(self.supplement_dir):
            try:
                os.mkdir(self.supplement_dir)
            except Exception as e:
                raise Exception("Output directory %s could not be created: %s" % (self.supplement_dir, str(e)))

        if self.is_phantomjs_installed():
            self.supplement_png_dir = os.path.join(self.supplement_dir, "png")
            if not os.path.exists(self.supplement_png_dir):
                try:
                    os.mkdir(os.path.join(self.supplement_dir, "png"))
                except Exception as e:
                    raise Exception("Output directory %s could not be created: %s" % (self.supplement_png_dir, str(e)))

            self.supplement_svg_dir = os.path.join(self.supplement_dir, "svg")
            if not os.path.exists(os.path.join(self.supplement_dir, "svg")):
                try:
                    os.mkdir(os.path.join(self.supplement_dir, "svg"))
                except Exception as e:
                    raise Exception("Output directory %s could not be created: %s" % (self.supplement_svg_dir, str(e)))

    def initialize_data(self):
        self.logger.info("Reading data counts.")
        for d in self.count_folders:
            sample = os.path.basename(d)
            self.samples.append(sample)
            self.update_dicts(sample, d)

        self.colors_dict = dict(
            zip(self.samples,
                self.large_colors_list[0:len(self.samples)]))

        self.logger.info("Finished reading data counts.")

        if len(self.count_folders) > 1:
            try:
                self.combine_counts()
                self.combine = True
            except Exception as e:
                err_msg = "Could not combine data counts: %s. Combined plots will not be generated" % str(e)
                self.logger.warn(err_msg)

    def update_dicts(self, sample, directory):
        self.logger.info("reading counts for sample: %s" % sample)
        fivepseq_out = FivePSeqOut(directory)

        self.data_summary_dict.update({sample: self.read_data_summary(fivepseq_out)})
        self.meta_count_start_dict.update({sample: self.read_meta_count_start(fivepseq_out)})
        self.meta_count_term_dict.update({sample: self.read_meta_count_term(fivepseq_out)})

        self.frame_count_term_dict.update({sample: self.read_frame_count_term(fivepseq_out)})
        self.frame_count_start_dict.update({sample: self.read_frame_count_start(fivepseq_out)})
        self.frame_stats_df_dict.update({sample: self.read_frame_stats_df(fivepseq_out)})
        self.amino_acid_df_dict.update({sample: self.read_amino_acid_df(fivepseq_out, full=False)})
        self.amino_acid_df_full_dict.update({sample: self.read_amino_acid_df(fivepseq_out, full=True)})
        self.codon_df_dict.update({sample: self.read_codon_df(fivepseq_out, basesort=False)})
        self.codon_basesorted_df_dict.update({sample: self.read_codon_df(fivepseq_out, basesort=True)})

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
            amino_acid_df_full = self.amino_acid_df_full_dict.get(key)
            codon_df = self.codon_df_dict.get(key)
            codon_basesorted_df = self.codon_basesorted_df_dict.get(key)
            if len(self.meta_count_start_combined) == 0:
                self.meta_count_start_combined = df_start.copy()
                self.meta_count_term_combined = df_term.copy()
                self.frame_count_START_combined = frame_start.copy()
                self.frame_count_TERM_combined = frame_term.copy()
                self.amino_acid_df_combined = amino_acid_df.copy()
                self.amino_acid_df_full_combined = amino_acid_df_full.copy()
                self.codon_df_combined = codon_df.copy()
                self.codon_basesorted_df_combined = codon_basesorted_df.copy()
            else:
                self.meta_count_start_combined.C += df_start.C
                self.meta_count_term_combined.C += df_term.C
                self.frame_count_START_combined.loc[:, ('F0', 'F1', 'F2')] += frame_start.loc[:, ('F0', 'F1', 'F2')]
                self.frame_count_TERM_combined.loc[:, ('F0', 'F1', 'F2')] += frame_term.loc[:, ('F0', 'F1', 'F2')]
                if amino_acid_df is not None:
                    self.amino_acid_df_combined += amino_acid_df
                if amino_acid_df_full is not None:
                    self.amino_acid_df_full_combined += amino_acid_df_full
                if codon_df is not None:
                    self.codon_df_combined += codon_df
                if codon_basesorted_df is not None:
                    self.codon_basesorted_df_combined += codon_basesorted_df

    def plot_multiple_samples(self):
        self.logger.info("Generating plots")

        self.make_single_sample_plots(self.title)

        if self.combine:
            self.make_combined_plots(self.title)

            bokeh_composite(self.title,
                            [self.p_scatter_start, self.p_scatter_term, None, None,
                             self.p_scatter_start_scaled, self.p_scatter_term_scaled, None, None,
                             self.p_triangle_term, self.p_triangle_start, None, None,
                             self.p_aa_heatmap, None, None, None,
                             self.p_aa_heatmap_scaled, None, None, None,
                             self.p_frame_barplots_term, None, None, None,
                             self.p_frame_barplots_start, None, None, None,
                             self.p_fft_plot_start, self.p_fft_plot_term, None, None,
                             self.p_scatter_start_combined, self.p_scatter_term_combined, None, None,
                             self.p_scatter_start_scaled_combined, self.p_scatter_term_scaled_combined, None, None,
                             self.p_triangle_term_combined, self.p_triangle_start_combined, None, None,
                             self.p_aa_heatmaps_combined, None, None, None],
                            os.path.join(self.args.o, self.title + ".html"), 4)
        else:
            bokeh_composite(self.title,
                            [self.p_scatter_start, self.p_scatter_term, None, None,
                             self.p_scatter_start_scaled, self.p_scatter_term_scaled, None, None,
                             self.p_triangle_term, self.p_triangle_start, None, None,
                             self.p_aa_heatmap, None, None, None,
                             self.p_aa_heatmap_scaled, None, None, None,
                             self.p_frame_barplots_term, None, None, None,
                             self.p_frame_barplots_start, None, None, None,
                             self.p_fft_plot_start, self.p_fft_plot_term, None, None],
                            os.path.join(self.args.o, self.title + ".html"), 4)

    def write_supplement(self):
        # codon pauses
        self.logger.info("Generating supplement plots: codon pauses")
        codon_title = self.title + "_codon_pauses"

        if self.combine:
            bokeh_composite(codon_title,
                            [bokeh_heatmap_grid(codon_title, self.codon_df_dict, scale=False),
                             bokeh_heatmap_grid(codon_title + "_scaled", self.codon_df_dict, scale=True),
                             bokeh_heatmap_grid(codon_title + "_combined", {"combined:": self.codon_df_combined},
                                                scale=False),
                             bokeh_heatmap_grid(codon_title + "_combined_scaled", {"combined:": self.codon_df_combined},
                                                scale=True),
                             bokeh_heatmap_grid(codon_title + "_basesorted", self.codon_basesorted_df_dict, scale=False),
                             bokeh_heatmap_grid(codon_title + "_basesorted_scaled", self.codon_basesorted_df_dict,
                                                scale=True),
                             bokeh_heatmap_grid(codon_title + "_basesorted_combined",
                                                {"basesored_combined": self.codon_basesorted_df_combined}, scale=False),
                             bokeh_heatmap_grid(codon_title + "_basesorted_combined_scaled",
                                                {"basesored_combined": self.codon_basesorted_df_combined}, scale=True)],
                            os.path.join(self.supplement_dir, codon_title + ".html"), 1)
        else:
            bokeh_composite(codon_title,
                            [bokeh_heatmap_grid(codon_title, self.codon_df_dict, scale=False),
                             bokeh_heatmap_grid(codon_title + "_scaled", self.codon_df_dict, scale=True),
                             bokeh_heatmap_grid(codon_title + "_basesorted", self.codon_basesorted_df_dict,
                                                scale=False),
                             bokeh_heatmap_grid(codon_title + "_basesorted_scaled", self.codon_basesorted_df_dict,
                                                scale=True)],
                            os.path.join(self.supplement_dir, codon_title + ".html"), 1)
        # amino acid scatter-plots
        self.logger.info("Generating supplement plots: amino acid scatter-plots")
        aa_scatterplots = []
        for aa in Codons.AMINO_ACID_TABLE.keys():
            self.logger.info("Plotting scatter for %s counts" % aa)

            aa_count_dict = {}

            for sample in self.samples:
                amino_acid_df_full = self.amino_acid_df_full_dict[sample]
                aa_df = pd.DataFrame(data={'D': map(int, amino_acid_df_full.columns), 'C': amino_acid_df_full.loc[aa, :]})
                aa_df = aa_df.reset_index(drop=True)
                aa_count_dict.update({sample: aa_df})
            aa_sp = bokeh_scatter_plot(aa, FivePSeqCounts.TERM, aa_count_dict, self.colors_dict, scale=True,
                               png_dir=self.supplement_png_dir, svg_dir=self.supplement_svg_dir)
            aa_scatterplots.append(aa_sp)

            if self.combine:
                aa_df = pd.DataFrame(data={'D': map(int, self.amino_acid_df_full_combined.columns),
                                           'C': self.amino_acid_df_full_combined.loc[aa, :]})
                aa_df = aa_df.reset_index(drop=True)
                aa_sp = bokeh_scatter_plot(aa + "_combined", FivePSeqCounts.TERM, {"combined": aa_df}, self.combined_color_dict,
                                   scale=True,
                                   png_dir=self.supplement_png_dir, svg_dir=self.supplement_svg_dir)

                aa_scatterplots.append(aa_sp)

        bokeh_composite(self.title + "amino_acid_scatterplots", aa_scatterplots,
                        os.path.join(self.supplement_dir, self.title + "amino_acid_scatterplots.html"), 2)

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

    def read_amino_acid_df(self, fivepseq_out, full=False):
        file = fivepseq_out.get_file_path(FivePSeqOut.AMINO_ACID_PAUSES_FILE)
        try:
            amino_acid_df = CountManager.read_amino_acid_df(file)

            if not full:
                if amino_acid_df.shape[1] > self.dist_for_amino_acid_heatmaps:
                    colrange = map(str, np.arange(-1 * self.dist_for_amino_acid_heatmaps, 0))
                    amino_acid_df = amino_acid_df.loc[:, colrange]

        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            amino_acid_df = None
        return amino_acid_df

    def read_codon_df(self, fivepseq_out, basesort=False):
        file = fivepseq_out.get_file_path(FivePSeqOut.CODON_PAUSES_FILE)
        try:
            codon_df = CountManager.read_amino_acid_df(file)
            if basesort:
                sorted_index = [""] * len(codon_df.index)
                for i in range(len(codon_df.index)):
                    ind = codon_df.index[i]
                    aa = ind.split("_")[0]
                    codon = ind.split("_")[1]
                    sorted_index[i] = codon + "_" + aa
                codon_df.index = sorted_index
                codon_df = codon_df.reindex(sorted(sorted_index))
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            codon_df = None
        return codon_df

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

        self.p_scatter_term_scaled = bokeh_scatter_plot(title + self.METACOUNTS_TERM_SCALED, FivePSeqCounts.TERM,
                                                        self.meta_count_term_dict,
                                                        self.colors_dict,
                                                        scale=True,
                                                        png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_scatter_start_scaled = bokeh_scatter_plot(title + self.METACOUNTS_START_SCALED, FivePSeqCounts.START,
                                                         self.meta_count_start_dict,
                                                         self.colors_dict,
                                                         scale=True,
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
        self.p_scatter_term_combined = bokeh_scatter_plot(title + self.METACOUNTS_TERM + "_combined",
                                                          FivePSeqCounts.TERM,
                                                          {"combined": self.meta_count_term_combined},
                                                          self.combined_color_dict,
                                                          png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_scatter_start_combined = bokeh_scatter_plot(title + self.METACOUNTS_START + "_combined",
                                                           FivePSeqCounts.START,
                                                           {self.COMBINED: self.meta_count_start_combined},
                                                           self.combined_color_dict,
                                                           png_dir=self.png_dir, svg_dir=self.svg_dir)

        # Combined plots
        self.p_scatter_term_scaled_combined = bokeh_scatter_plot(title + self.METACOUNTS_TERM_SCALED + "_combined",
                                                                 FivePSeqCounts.TERM,
                                                                 {"combined": self.meta_count_term_combined},
                                                                 self.combined_color_dict,
                                                                 scale=True,
                                                                 png_dir=self.png_dir, svg_dir=self.svg_dir)

        self.p_scatter_start_scaled_combined = bokeh_scatter_plot(title + self.METACOUNTS_START_SCALED + "_combined",
                                                                  FivePSeqCounts.START,
                                                                  {self.COMBINED: self.meta_count_start_combined},
                                                                  self.combined_color_dict,
                                                                  scale=True,
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
