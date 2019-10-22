import glob
import logging
import os

import colorlover as cl
import numpy as np
import pandas as pd
from bokeh.io import export_svgs
from bokeh.models import Div
from bokeh.plotting import figure

from fivepseq import config
from fivepseq.logic.structures.codons import Codons
from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from fivepseq.util.writers import FivePSeqOut
from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_heatmap_grid, bokeh_frame_barplots, \
    bokeh_composite, bokeh_fft_plot, bokeh_tabbed_scatter_plot, bokeh_tabbed_triangle_plot, \
    bokeh_tabbed_frame_barplots, bokeh_tabbed_heatmap_grid, bokeh_tabbed_fft_plot
from fivepseq.viz.header_html import write_fivepseq_header


class VizPipeline:
    FILTER_TOP_POPULATED = "populated"
    FILTER_CANONICAL_TRANSCRIPTS = "canonical"

    DATA_TYPE_META_COUNTS = 0
    DATA_TYPE_FRAME_COUNTS = 1
    DATA_TYPE_FFT_SIGNALS = 2
    DATA_TYPE_AMINO_ACID_DF = 3
    DATA_TYPE_FRAME_STATS = 4
    DATA_TYPE_LIB_SIZE = 5

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
    sample_folders_dict = {}
    title = "fivepseq_plot_canvas"
    main_dir = None
    png_dir = None
    svg_dir = None
    supplement_dir = "supplement"
    supplement_png_dir = None
    supplement_svg_dir = None

    geneset_dir = "genesets"
    geneset_png_dir = None
    geneset_svg_dir = None

    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

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
    loci_meta_counts_dict_ALL = {}
    loci_meta_counts_dict_3UTR = {}
    loci_meta_counts_dict_5UTR = {}
    loci_meta_counts_dict_CDS = {}
    data_summary_dict = {}
    lib_size_dict = {}
    transcript_index = None

    # the distance to be used for plotting amino-acid heatmaps
    dist_for_amino_acid_heatmaps = 20

    gs_transcriptInd_dict = None
    gs_shortGS_dict = {}
    shortGS_gs_dict = {}

    # large_colors_list = (cl.to_numeric(cl.scales['8']['qual']['Paired'][0:6]) +
    #                     cl.to_numeric(cl.scales['8']['qual']['Set1'][3:5]))
    large_colors_list = ([cl.scales['8']['qual']['Set1'][i] for i in [1, 2, 3, 4, 0]] +
                         [cl.scales['8']['qual']['Paired'][i] for i in [0, 2, 4]])
    #                     cl.to_numeric(cl.scales['5']['qual']['Set3']))
    gs_colors_list = ("#AA4488", "#4477AA", "#AAAA44", "#AA7744", "#AA4455", "#44AAAA",
                      "#771155", "#114477", "#777711", "#774411", "#771122", "#117777")

    colors_dict = None
    # combined_color_dict = {COMBINED: cl.to_numeric(cl.scales['9']['qual']['Set3'])[3]}
    COMBINED = "combined"

    combine_sum_color = cl.scales['3']['div']['RdGy'][0]
    combine_weighted_color = cl.scales['3']['div']['BrBG'][0]

    combined_color_dict = {COMBINED: combine_sum_color}
    combined_weighted_color_dict = {COMBINED: combine_weighted_color}

    phantomjs_installed = None

    fivepseq_counts_container_dict = {}

    combine = False  # do not combine plots until data counts are successfully combined

    figure_list = None
    figure_list_combined = None

    def __init__(self, args, count_folders=None, gs_transcriptInd_dict=None):
        """
        Initialize vizualization pipeline with arguments contained in args.
        If count_folders are provided explicitely, those will be used instead of sd and md arguments.

        :param args:
        :param count_folders:
        """

        self.gs_transcriptInd_dict = gs_transcriptInd_dict
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
            self.plot_main()
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

        try:
            if self.gs_transcriptInd_dict is not None:
                self.write_genesets()
        except Exception as e:
            err_msg = "Exception while writing genesets: %s" % str(e)
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

            elif hasattr(self.args, 'md') and self.args.md is not None:
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

        self.main_dir = os.path.join(self.args.o, "main")
        if not os.path.exists(self.main_dir):
            try:
                os.mkdir(self.main_dir)
            except Exception as e:
                raise Exception("Output directory %s could not be created: %s" % (self.args.o, str(e)))

        if self.is_phantomjs_installed():
            self.png_dir = os.path.join(self.main_dir, "png")
            if not os.path.exists(self.png_dir):
                try:
                    os.mkdir(self.png_dir)
                except Exception as e:
                    raise Exception("Output directory %s could not be created: %s" % (self.png_dir, str(e)))

            self.svg_dir = os.path.join(self.main_dir, "svg")
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

        if self.gs_transcriptInd_dict is not None:
            self.geneset_dir = os.path.join(self.args.o, "genesets")
            if not os.path.exists(self.geneset_dir):
                try:
                    os.mkdir(self.geneset_dir)
                except Exception as e:
                    raise Exception("Output directory %s could not be created: %s" % (self.geneset_dir, str(e)))

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

            if self.gs_transcriptInd_dict is not None:
                self.geneset_png_dir = os.path.join(self.geneset_dir, "png")
                if not os.path.exists(self.geneset_png_dir):
                    try:
                        os.mkdir(os.path.join(self.geneset_dir, "png"))
                    except Exception as e:
                        raise Exception("Output directory %s could not be created: %s" % (self.geneset_png_dir, str(e)))

                self.geneset_svg_dir = os.path.join(self.geneset_dir, "svg")
                if not os.path.exists(os.path.join(self.geneset_dir, "svg")):
                    try:
                        os.mkdir(os.path.join(self.geneset_dir, "svg"))
                    except Exception as e:
                        raise Exception("Output directory %s could not be created: %s" % (self.geneset_svg_dir, str(e)))

    def initialize_data(self):
        self.logger.info("Reading data counts.")
        for d in self.count_folders:
            sample = os.path.basename(os.path.dirname(d))
            self.samples.append(sample)
            self.sample_folders_dict.update({sample: os.path.dirname(d)})
            self.update_dicts(sample, d)

        self.colors_dict = dict(
            zip(self.samples,
                self.large_colors_list[0:len(self.samples)]))

        self.logger.info("Finished reading data counts.")

        if len(self.count_folders) > 1:
            try:
                self.combine = True
            except Exception as e:
                err_msg = "Could not combine data counts: %s. Combined plots will not be generated" % str(e)
                self.logger.warn(err_msg)

    def update_dicts(self, sample, directory):
        self.logger.info("reading counts for sample: %s" % sample)
        fivepseq_out = FivePSeqOut(directory)

        data_summary = self.read_data_summary(fivepseq_out)
        self.data_summary_dict.update({sample: data_summary})
        self.lib_size_dict.update({sample: data_summary.iloc[0, 0]})
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

        self.loci_meta_counts_dict_ALL.update({sample: self.read_loci_meta_counts(
            fivepseq_out, file=fivepseq_out.get_file_path(
                FivePSeqOut.LOCI_PAUSES_FILE_PREFIX + FivePSeqCounts.READ_LOCATIONS_ALL + ".txt"))})
        self.loci_meta_counts_dict_3UTR.update({sample: self.read_loci_meta_counts(
            fivepseq_out, file=fivepseq_out.get_file_path(
                FivePSeqOut.LOCI_PAUSES_FILE_PREFIX + FivePSeqCounts.READ_LOCATIONS_3UTR + ".txt"))})
        self.loci_meta_counts_dict_5UTR.update({sample: self.read_loci_meta_counts(
            fivepseq_out, file=fivepseq_out.get_file_path(
                FivePSeqOut.LOCI_PAUSES_FILE_PREFIX + FivePSeqCounts.READ_LOCATIONS_5UTR + ".txt"))})
        self.loci_meta_counts_dict_CDS.update({sample: self.read_loci_meta_counts(
            fivepseq_out, file=fivepseq_out.get_file_path(
                FivePSeqOut.LOCI_PAUSES_FILE_PREFIX + FivePSeqCounts.READ_LOCATIONS_CDS + ".txt"))})

        # TODO remove transcript filter
        if hasattr(self.args, "tf") and self.args.tf is not None:
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

    def plot_main(self):
        self.logger.info("Generating plots")
        # plots header page
        # bokeh_composite(self.title + "_main",
        #                header,
        #                os.path.join(self.args.o, self.title + ".html"), 1)
        write_fivepseq_header(self)

        self.figure_list = []
        self.figure_list += [self.get_scatter_plot(region=FivePSeqCounts.START),
                             self.get_scatter_plot(region=FivePSeqCounts.TERM)]
        self.figure_list += [self.get_scatter_plot(region=FivePSeqCounts.START, scale=True),
                             self.get_scatter_plot(region=FivePSeqCounts.TERM, scale=True)]
        self.figure_list += [self.get_frame_barplot(), None]
        self.figure_list += [self.get_triangle_plot(), None]  # self.p_triangle_start]
        self.figure_list += [self.get_fft_plot(region=FivePSeqCounts.START),
                             self.get_fft_plot(region=FivePSeqCounts.TERM)]
        self.figure_list += [self.get_heatmap_plot(), None]
        self.figure_list += [self.get_heatmap_plot(scale=True), None]

        if self.loci_meta_counts_dict_ALL is not None:
            self.figure_list += [self.get_scatter_plot(plot_name="ALL_metacounts_relative_to", region="loci",
                                                       count_dict=self.loci_meta_counts_dict_ALL, scale=True)]
        if self.loci_meta_counts_dict_3UTR is not None:
            self.figure_list += [self.get_scatter_plot(plot_name="ALL_metacounts_relative_to", region="loci",
                                                       count_dict=self.loci_meta_counts_dict_3UTR, scale=True)]
        if self.loci_meta_counts_dict_5UTR is not None:
            self.figure_list += [self.get_scatter_plot(plot_name="5UTR_metacounts_relative_to", region="loci",
                                                       count_dict=self.loci_meta_counts_dict_5UTR, scale=True)]
        if self.loci_meta_counts_dict_CDS is not None:
            self.figure_list += [self.get_scatter_plot(plot_name="CDS_metacounts_relative_to", region="loci",
                                                       count_dict=self.loci_meta_counts_dict_CDS, scale=True)]

        bokeh_composite(self.title + "_main",
                        self.figure_list,
                        os.path.join(self.main_dir, self.title + "_main.html"), 2)

        if self.combine:
            self.figure_list_combined = [self.get_scatter_plot(region=FivePSeqCounts.START,
                                                               combine_sum=True,
                                                               combine_color=self.combine_sum_color),
                                         self.get_scatter_plot(region=FivePSeqCounts.TERM,
                                                               combine_sum=True,
                                                               combine_color=self.combine_sum_color)]
            self.figure_list_combined += [self.get_scatter_plot(region=FivePSeqCounts.START,
                                                                combine_weighted=True,
                                                                combine_color=self.combine_weighted_color),
                                          self.get_scatter_plot(region=FivePSeqCounts.TERM,
                                                                combine_weighted=True,
                                                                combine_color=self.combine_weighted_color)]
            self.figure_list_combined += [
                self.get_scatter_plot(region=FivePSeqCounts.START, combine_sum=True,
                                      combine_color=self.combine_sum_color, scale=True),
                self.get_scatter_plot(region=FivePSeqCounts.TERM, combine_sum=True,
                                      combine_color=self.combine_sum_color, scale=True)]
            self.figure_list_combined += [
                self.get_scatter_plot(region=FivePSeqCounts.START, combine_weighted=True,
                                      combine_color=self.combine_weighted_color, scale=True),
                self.get_scatter_plot(region=FivePSeqCounts.TERM, combine_weighted=True,
                                      combine_color=self.combine_weighted_color, scale=True)]
            self.figure_list_combined += [
                self.get_triangle_plot(combine_sum=True, combine_color=self.combine_sum_color), None]
            self.figure_list_combined += [
                self.get_triangle_plot(combine_weighted=True, combine_color=self.combine_weighted_color), None]
            self.figure_list_combined += [self.get_heatmap_plot(scale=True, combine_sum=True), None]
            self.figure_list_combined += [self.get_heatmap_plot(scale=True, combine_weighted=True), None]

            bokeh_composite(self.title + "_" + self.COMBINED,
                            self.figure_list_combined,
                            os.path.join(self.main_dir, self.title + "_" + self.COMBINED + ".html"), 2)

    GS_COUNTER = 1

    def generate_shortGS(self):
        shortGS = "GS" + str(self.GS_COUNTER)
        self.GS_COUNTER += 1
        return shortGS

    def write_genesets(self):

        self.logger.info("Generating geneset-specific plots")

        for gs in self.gs_transcriptInd_dict.keys():
            shortGS = self.generate_shortGS()
            self.gs_shortGS_dict.update({gs: shortGS})
            self.shortGS_gs_dict.update({shortGS: gs})

        # tabs = genesets, overlay = samples
        geneset_title = self.title + "_samples_per_geneset"

        tabbed_plots = self.gs_tabbed_plots(geneset_title)
        bokeh_composite(geneset_title, tabbed_plots, os.path.join(self.geneset_dir, geneset_title + ".html"), 1)

        # tabs = samples, overlay = genesets
        sample_geneset_title = self.title + "_genesets_per_sample"

        sample_gs_tabbed_plots = self.sample_geneset_tabbed_plots(sample_geneset_title)
        bokeh_composite(sample_geneset_title, sample_gs_tabbed_plots,
                        os.path.join(self.geneset_dir, sample_geneset_title + ".html"), 1)

        self.logger.info("Wrote geneset-specific plots")

    def write_supplement(self):
        # codon pauses
        self.logger.info("Generating supplement plots: codon pauses")
        codon_title = self.title + "_codon_relative_counts"

        if self.combine:
            bokeh_composite(codon_title,
                            [self.get_heatmap_plot(codon_title, self.codon_df_dict, scale=False),
                             self.get_heatmap_plot(codon_title, self.codon_df_dict, scale=True),
                             self.get_heatmap_plot(codon_title, self.codon_df_dict, scale=True, combine_weighted=True),

                             self.get_heatmap_plot(codon_title + "_basesorted", self.codon_basesorted_df_dict,
                                                   scale=False),
                             self.get_heatmap_plot(codon_title + "_basesorted", self.codon_basesorted_df_dict,
                                                   scale=True),
                             self.get_heatmap_plot(codon_title + "_basesorted", self.codon_basesorted_df_dict,
                                                   scale=True, combine_weighted=True),
                             ],
                            os.path.join(self.supplement_dir, codon_title + ".html"), 1)
        else:
            bokeh_composite(codon_title,
                            [self.get_heatmap_plot(codon_title, self.codon_df_dict, scale=False),
                             self.get_heatmap_plot(codon_title, self.codon_df_dict, scale=True),
                             self.get_heatmap_plot(codon_title + "_basesorted", self.codon_basesorted_df_dict,
                                                   scale=False),
                             bokeh_heatmap_grid(codon_title + "_basesorted", self.codon_basesorted_df_dict,
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
                aa_df = pd.DataFrame(
                    data={'D': map(int, amino_acid_df_full.columns), 'C': amino_acid_df_full.loc[aa, :]})
                aa_df = aa_df.reset_index(drop=True)
                aa_count_dict.update({sample: aa_df})

            aa_scatterplots.append(self.get_scatter_plot(plot_name=aa, region="codon", count_dict=aa_count_dict,
                                                         scale=True,
                                                         png_dir=self.supplement_png_dir,
                                                         svg_dir=self.supplement_svg_dir))

            if self.combine:
                aa_scatterplots.append(self.get_scatter_plot(plot_name=aa, region="codon", count_dict=aa_count_dict,
                                                             scale=True,
                                                             combine_weighted=True,
                                                             combine_color=self.combine_weighted_color,
                                                             png_dir=self.supplement_png_dir,
                                                             svg_dir=self.supplement_svg_dir))

        bokeh_composite(self.title + "_amino_acid_scatterplots", aa_scatterplots,
                        os.path.join(self.supplement_dir, self.title + "_amino_acid_scatterplots.html"), 2)

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
            logging.getLogger(config.FIVEPSEQ_LOGGER).warning(
                "It seems like phantomjs is not installed no your system. "
                "Files may not be exported in svg and png formats, while html will still be available for viewing."
                "To install phantomjs, run 'conda install phantomjs selenium pillow'")
        return self.phantomjs_installed

    def read_meta_count_term(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.META_COUNT_TERM_FILE)

        try:
            meta_count_term = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            meta_count_term = None
        return meta_count_term

    def read_meta_count_start(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.META_COUNT_START_FILE)

        try:
            meta_count_start = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            meta_count_start = None
        return meta_count_start

    def read_amino_acid_df(self, fivepseq_out=None, file=None, full=False):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
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

    def read_codon_df(self, fivepseq_out=None, file=None, basesort=False):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
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

    def read_frame_count_term(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.FRAME_COUNTS_TERM_FILE)

        try:
            frame_count_term = CountManager.read_frame_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            frame_count_term = None
        return frame_count_term

    def read_frame_count_start(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.FRAME_COUNTS_START_FILE)
        try:
            frame_count_start = CountManager.read_frame_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            frame_count_start = None
        return frame_count_start

    def read_frame_stats_df(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.FRAME_STATS_DF_FILE)

        if os.path.exists(file):
            frame_stats_df = pd.read_csv(file, sep="\t", header=0, index_col=0)
        else:
            frame_stats_df = None
        return frame_stats_df

    def read_count_vector_list_term(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.COUNT_TERM_FILE)

        try:
            count_vector_list_term = CountManager.read_counts_as_list(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            count_vector_list_term = None
        return count_vector_list_term

    def read_count_vector_list_start(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.COUNT_START_FILE)

        try:
            count_vector_list_start = CountManager.read_counts_as_list(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            count_vector_list_start = None
        return count_vector_list_start

    def read_loci_meta_counts(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.LOCI_PAUSES_FILE)

        loci_meta_counts = None
        if os.path.exists(file):
            self.logger.info("Loci count file found:" + str(file))
            loci_meta_counts = CountManager.read_meta_counts(file)
        return loci_meta_counts

    def read_fft_signal_start(self, fivepseq_out, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.FFT_SIGNALS_START)
        try:
            fft_signals_start = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            fft_signals_start = None
        return fft_signals_start

    def read_fft_signal_term(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.FFT_SIGNALS_TERM)
        try:
            fft_signals_term = CountManager.read_meta_counts(file)
        except:
            self.logger.warn("The file %s not found, plots for this will be skipped." % str(file))
            fft_signals_term = None
        return fft_signals_term

    def read_data_summary(self, fivepseq_out=None, file=None):
        if file is None:
            if fivepseq_out is None:
                raise Exception("Insufficient arguments")
            file = fivepseq_out.get_file_path(FivePSeqOut.DATA_SUMMARY_FILE)
        try:
            data_summary = pd.read_csv(file, sep="\t", header=None, index_col=0)
        except:
            self.logger.warn("The file %s was not found, table will not be generated" % str(file))
            data_summary = None
        return data_summary

    def generate_gs_sample_dict(self, filename, data_type):
        gs_sample_dict = {}
        for gs in self.gs_transcriptInd_dict.keys():
            shortGS = self.gs_shortGS_dict[gs]
            gs_sample_dict.update({shortGS: {}})
            for sample in self.samples:
                if data_type == self.DATA_TYPE_META_COUNTS:
                    gs_sample_dict[shortGS].update({sample: self.read_meta_count_term(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_FRAME_COUNTS:
                    gs_sample_dict[shortGS].update({sample: self.read_frame_count_term(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_FRAME_STATS:
                    gs_sample_dict[shortGS].update({sample: self.read_frame_stats_df(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_AMINO_ACID_DF:
                    gs_sample_dict[shortGS].update({sample: self.read_amino_acid_df(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_FFT_SIGNALS:
                    gs_sample_dict[shortGS].update({sample: self.read_fft_signal_term(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_LIB_SIZE:
                    gs_sample_dict[shortGS].update({sample: self.lib_size_dict.get(sample)})

        return gs_sample_dict

    def generate_sample_gs_dict(self, filename, data_type):
        sample_gs_dict = {}
        for sample in self.samples:
            sample_gs_dict.update({sample: {}})
            for gs in self.gs_transcriptInd_dict.keys():
                shortGS = self.gs_shortGS_dict[gs]
                if data_type == self.DATA_TYPE_META_COUNTS:
                    sample_gs_dict[sample].update({shortGS: self.read_meta_count_term(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_FRAME_COUNTS:
                    sample_gs_dict[sample].update({shortGS: self.read_frame_count_term(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_FRAME_STATS:
                    sample_gs_dict[sample].update({shortGS: self.read_frame_stats_df(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_AMINO_ACID_DF:
                    sample_gs_dict[sample].update({shortGS: self.read_amino_acid_df(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_FFT_SIGNALS:
                    sample_gs_dict[sample].update({shortGS: self.read_fft_signal_term(
                        file=os.path.join(self.sample_folders_dict[sample], gs, filename))})
                elif data_type == self.DATA_TYPE_LIB_SIZE:
                    sample_gs_dict[sample].update({shortGS: self.lib_size_dict.get(sample)})
        return sample_gs_dict

    def gs_tabbed_plots(self, title):

        plots = []

        # gs_sample_count_series_dict_dict
        gs_meta_count_start_dict = self.generate_gs_sample_dict(FivePSeqOut.META_COUNT_START_FILE,
                                                                self.DATA_TYPE_META_COUNTS)
        gs_meta_count_term_dict = self.generate_gs_sample_dict(FivePSeqOut.META_COUNT_TERM_FILE,
                                                               self.DATA_TYPE_META_COUNTS)
        gs_frame_count_term_dict = self.generate_gs_sample_dict(FivePSeqOut.FRAME_COUNTS_TERM_FILE,
                                                                self.DATA_TYPE_FRAME_COUNTS)
        gs_frame_stats_dict = self.generate_gs_sample_dict(FivePSeqOut.FRAME_STATS_DF_FILE, self.DATA_TYPE_FRAME_STATS)
        gs_amino_acid_df_dict = self.generate_gs_sample_dict(FivePSeqOut.AMINO_ACID_PAUSES_FILE,
                                                             self.DATA_TYPE_AMINO_ACID_DF)
        gs_fft_start_dict = self.generate_gs_sample_dict(FivePSeqOut.FFT_SIGNALS_START, self.DATA_TYPE_FFT_SIGNALS)
        gs_fft_term_dict = self.generate_gs_sample_dict(FivePSeqOut.FFT_SIGNALS_TERM, self.DATA_TYPE_FFT_SIGNALS)
        gs_lib_size_dict = self.generate_gs_sample_dict(None, self.DATA_TYPE_LIB_SIZE)

        plots.append(self.get_gs_mapping_div())
        plots.append(self.get_tabbed_scatter_plot(group_count_dict=gs_meta_count_start_dict,
                                                  region=FivePSeqCounts.START,
                                                  scale=True, lib_size_dict_dict=gs_lib_size_dict,
                                                  png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_scatter_plot(group_count_dict=gs_meta_count_term_dict,
                                                  region=FivePSeqCounts.TERM,
                                                  scale=True, lib_size_dict_dict=gs_lib_size_dict,
                                                  png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_frame_barplots(group_count_dict=gs_frame_count_term_dict,
                                                    group_frame_stats_df_dict=gs_frame_stats_dict,
                                                    png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_triangle_plot(group_count_dict=gs_frame_count_term_dict,
                                                   png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_fft_plot(group_count_dict=gs_fft_start_dict,
                                              region=FivePSeqCounts.START,
                                              png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_fft_plot(group_count_dict=gs_fft_start_dict,
                                              region=FivePSeqCounts.TERM,
                                              png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_heatmap_plot(group_count_dict=gs_amino_acid_df_dict, scale=True,
                                                  png_dir=self.geneset_png_dir,
                                                  svg_dir=self.geneset_svg_dir))

        # combined dicts
        if self.combine:
            gs_meta_counts_start_combined = {}
            gs_meta_counts_term_combined = {}
            gs_frame_counts_term_combined = {}
            gs_amino_acid_counts_combined = {}
            gs_lib_size_combined = {}

            for gs in self.shortGS_gs_dict.keys():
                gs_meta_counts_start_combined.update(
                    {gs: {self.COMBINED: CountManager.combine_count_series(gs_meta_count_start_dict[gs])}})

                gs_meta_counts_term_combined.update(
                    {gs: {self.COMBINED: CountManager.combine_count_series(gs_meta_count_term_dict[gs])}})

                gs_frame_counts_term_combined.update(
                    # {gs: {self.COMBINED: self.combine_frame_counts(gs_frame_count_term_dict[gs])}})
                    {gs: {self.COMBINED: CountManager.combine_frame_counts(gs_frame_count_term_dict[gs],
                                                                           gs_lib_size_dict[gs])}})

                gs_amino_acid_counts_combined.update(
                    {gs: {self.COMBINED: CountManager.combine_amino_acid_dfs(gs_amino_acid_df_dict[gs],
                                                                             gs_lib_size_dict[gs])}})

                gs_lib_size_combined.update({gs: {self.COMBINED: sum(self.lib_size_dict.values())}})

            plots.append(self.get_tabbed_scatter_plot(group_count_dict=gs_meta_count_start_dict,
                                                      region=FivePSeqCounts.START,
                                                      scale=True, lib_size_dict_dict=gs_lib_size_dict,
                                                      combine_weighted=True, combine_color=self.combine_weighted_color,
                                                      png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

            plots.append(self.get_tabbed_scatter_plot(group_count_dict=gs_meta_count_term_dict,
                                                      region=FivePSeqCounts.TERM,
                                                      scale=True, lib_size_dict_dict=gs_lib_size_dict,
                                                      combine_weighted=True, combine_color=self.combine_weighted_color,
                                                      png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

            plots.append(self.get_tabbed_triangle_plot(group_count_dict=gs_frame_count_term_dict,
                                                       combine_weighted=True, combine_color=self.combine_weighted_color,
                                                       lib_size_dict_dict=gs_lib_size_dict,
                                                       png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

            plots.append(
                self.get_tabbed_heatmap_plot(group_count_dict=gs_amino_acid_counts_combined, scale=True,
                                             png_dir=self.geneset_png_dir,
                                             svg_dir=self.geneset_svg_dir))

        return plots

    def get_gs_mapping_div(self):
        table = "<table>" \
                "<tr>" \
                "<th>Gene-set actual name</th>" \
                "<th>Gene-set short name</th>" \
                "</tr>"

        for gs in self.gs_shortGS_dict.keys():
            table += "<tr><td>%s</td><td>%s</td></tr>" \
                     % (gs, self.gs_shortGS_dict[gs])
        table += "</table>"

        div_text = """<hr><div>
            <h4>Gene-set name mapping</h4>
            %s
            </div>""" % table
        div = Div(text=div_text)
        return div

    def sample_geneset_tabbed_plots(self, title):
        plots = []

        # gs_sample_count_series_dict_dict
        sample_gs_meta_count_start_dict = self.generate_sample_gs_dict(FivePSeqOut.META_COUNT_START_FILE,
                                                                       self.DATA_TYPE_META_COUNTS)
        sample_gs_meta_count_term_dict = self.generate_sample_gs_dict(FivePSeqOut.META_COUNT_TERM_FILE,
                                                                      self.DATA_TYPE_META_COUNTS)
        sample_gs_frame_count_term_dict = self.generate_sample_gs_dict(FivePSeqOut.FRAME_COUNTS_TERM_FILE,
                                                                       self.DATA_TYPE_FRAME_COUNTS)
        sample_gs_frame_stats_dict = self.generate_sample_gs_dict(FivePSeqOut.FRAME_STATS_DF_FILE,
                                                                  self.DATA_TYPE_FRAME_STATS)
        sample_gs_amino_acid_df_dict = self.generate_sample_gs_dict(FivePSeqOut.AMINO_ACID_PAUSES_FILE,
                                                                    self.DATA_TYPE_AMINO_ACID_DF)
        sample_gs_fft_start_dict = self.generate_sample_gs_dict(FivePSeqOut.FFT_SIGNALS_START,
                                                                self.DATA_TYPE_FFT_SIGNALS)
        sample_gs_fft_term_dict = self.generate_sample_gs_dict(FivePSeqOut.FFT_SIGNALS_TERM, self.DATA_TYPE_FFT_SIGNALS)

        sample_gs_lib_size_dict = self.generate_sample_gs_dict(None, self.DATA_TYPE_LIB_SIZE)

        gs_colors_dict = dict(
            zip(self.shortGS_gs_dict.keys(),
                self.gs_colors_list[0:len(self.shortGS_gs_dict.keys())]))

        plots.append(self.get_gs_mapping_div())
        plots.append(self.get_tabbed_scatter_plot(group_count_dict=sample_gs_meta_count_start_dict,
                                                  region=FivePSeqCounts.START, color_dict=gs_colors_dict,
                                                  scale=True, lib_size_dict_dict=sample_gs_lib_size_dict,
                                                  png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_scatter_plot(group_count_dict=sample_gs_meta_count_term_dict,
                                                  region=FivePSeqCounts.TERM, color_dict=gs_colors_dict,
                                                  scale=True, lib_size_dict_dict=sample_gs_lib_size_dict,
                                                  png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_frame_barplots(group_count_dict=sample_gs_frame_count_term_dict,
                                                    group_frame_stats_df_dict=sample_gs_frame_stats_dict,
                                                    color_dict=gs_colors_dict,
                                                    png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_triangle_plot(group_count_dict=sample_gs_frame_count_term_dict,
                                                   color_dict=gs_colors_dict,
                                                   lib_size_dict_dict=sample_gs_lib_size_dict,
                                                   png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_fft_plot(region=FivePSeqCounts.START,
                                              group_count_dict=sample_gs_fft_start_dict,
                                              color_dict=gs_colors_dict,
                                              png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_fft_plot(region=FivePSeqCounts.TERM,
                                              group_count_dict=sample_gs_fft_start_dict,
                                              color_dict=gs_colors_dict,
                                              png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

        plots.append(self.get_tabbed_heatmap_plot(group_count_dict=sample_gs_amino_acid_df_dict, scale=True,
                                                  png_dir=self.geneset_png_dir,
                                                  svg_dir=self.geneset_svg_dir))

        # combined plots
        if self.combine:
            # combined dicts
            sample_gs_meta_counts_start_combined = {}
            sample_gs_meta_counts_term_combined = {}
            sample_gs_frame_counts_term_combined = {}
            sample_gs_amino_acid_counts_combined = {}
            sample_gs_lib_size_combined = {}

            for gs in self.shortGS_gs_dict.keys():
                sample_temp_dict = {}
                for sample in self.samples:
                    sample_temp_dict.update({sample: sample_gs_meta_count_start_dict[sample][gs]})
                sample_gs_meta_counts_start_combined.update(
                    {gs: CountManager.combine_count_series(sample_temp_dict, self.lib_size_dict)})

                sample_temp_dict = {}
                for sample in self.samples:
                    sample_temp_dict.update({sample: sample_gs_meta_count_term_dict[sample][gs]})
                sample_gs_meta_counts_term_combined.update(
                    {gs: CountManager.combine_count_series(sample_temp_dict, self.lib_size_dict)})

                sample_temp_dict = {}
                for sample in self.samples:
                    sample_temp_dict.update({sample: sample_gs_frame_count_term_dict[sample][gs]})
                sample_gs_frame_counts_term_combined.update(
                    {gs: CountManager.combine_frame_counts(sample_temp_dict, self.lib_size_dict)})

                sample_temp_dict = {}
                for sample in self.samples:
                    sample_temp_dict.update({sample: sample_gs_amino_acid_df_dict[sample][gs]})
                sample_gs_amino_acid_counts_combined.update(
                    {gs: CountManager.combine_amino_acid_dfs(sample_temp_dict, self.lib_size_dict)})

                sample_gs_lib_size_combined.update({gs: sum(self.lib_size_dict.values())})

            plots.append(self.get_scatter_plot(region=FivePSeqCounts.START,
                                               count_dict=sample_gs_meta_counts_start_combined,
                                               color_dict=gs_colors_dict,
                                               scale=True, lib_size_dict=sample_gs_lib_size_combined,
                                               png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

            plots.append(self.get_scatter_plot(region=FivePSeqCounts.TERM,
                                               count_dict=sample_gs_meta_counts_term_combined,
                                               color_dict=gs_colors_dict,
                                               scale=True, lib_size_dict=sample_gs_lib_size_combined,
                                               png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))
            plots.append(self.get_triangle_plot(count_dict=sample_gs_frame_counts_term_combined,
                                                color_dict=gs_colors_dict,
                                                lib_size_dict=sample_gs_lib_size_combined,
                                                png_dir=self.geneset_png_dir, svg_dir=self.geneset_svg_dir))

            plots.append(self.get_heatmap_plot(count_dict=sample_gs_amino_acid_counts_combined, scale=True,
                                               png_dir=self.geneset_png_dir,
                                               svg_dir=self.geneset_svg_dir))

        return plots

    def get_scatter_plot(self, plot_name="metagene_counts", region="", count_dict=None,
                         color_dict=colors_dict,
                         scale=False, lib_size_dict=None,
                         combine_sum=False, combine_weighted=False, combine_color=None,
                         png_dir=png_dir, svg_dir=svg_dir):

        if count_dict is None:
            if region == FivePSeqCounts.TERM:
                count_dict = self.meta_count_term_dict
            else:
                count_dict = self.meta_count_start_dict

        if color_dict is None:
            color_dict = self.colors_dict

        if lib_size_dict is None:
            lib_size_dict = self.lib_size_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name + "-" + region

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        if scale:
            title += "_" + "scaled"
        else:
            title += "_" + "raw"

        p = bokeh_scatter_plot(title=title,
                               region=region,
                               count_series_dict=count_dict,
                               color_dict=color_dict,
                               scale=scale, lib_size_dict=lib_size_dict,
                               combine_sum=combine_sum,
                               combine_weighted=combine_weighted,
                               combine_color=combine_color,
                               png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_tabbed_scatter_plot(self, group_count_dict,
                                plot_name="metagene_counts", region="",
                                color_dict=colors_dict,
                                scale=False, lib_size_dict_dict=None,
                                combine_sum=False, combine_weighted=False, combine_color=None,
                                png_dir=png_dir, svg_dir=svg_dir):

        if color_dict is None:
            color_dict = self.colors_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name + "-" + region
        if scale:
            title += "_" + "scaled"
        else:
            title += "_" + "raw"

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_tabbed_scatter_plot(title=title,
                                      region=region,
                                      group_count_series_dict_dict=group_count_dict,
                                      color_dict=color_dict,
                                      scale=scale, lib_size_dict_dict=lib_size_dict_dict,
                                      combine_sum=combine_sum,
                                      combine_weighted=combine_weighted,
                                      combine_color=combine_color,
                                      png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_triangle_plot(self, plot_name="gene_frame_preferences",
                          count_dict=None,
                          color_dict=None, lib_size_dict=None,
                          combine_sum=False, combine_weighted=False, combine_color=None,
                          png_dir=png_dir, svg_dir=svg_dir):

        if count_dict is None:
            count_dict = self.frame_count_start_dict

        if color_dict is None:
            color_dict = self.colors_dict

        if lib_size_dict is None:
            lib_size_dict = self.lib_size_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_triangle_plot(title=title,
                                frame_df_dict=count_dict,
                                lib_size_dict=lib_size_dict,
                                color_dict=color_dict,
                                combine_sum=combine_sum,
                                combine_weighted=combine_weighted,
                                combine_color=combine_color,
                                png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_tabbed_triangle_plot(self, group_count_dict,
                                 plot_name="gene_frame_preferences",
                                 color_dict=None, lib_size_dict_dict=None,
                                 combine_sum=False, combine_weighted=False, combine_color=None,
                                 png_dir=png_dir, svg_dir=svg_dir):

        if color_dict is None:
            color_dict = self.colors_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_tabbed_triangle_plot(title=title,
                                       group_frame_df_dict_dict=group_count_dict,
                                       color_dict=color_dict,
                                       lib_size_dict_dict=lib_size_dict_dict,
                                       combine_sum=combine_sum,
                                       combine_weighted=combine_weighted,
                                       combine_color=combine_color,
                                       png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_heatmap_plot(self, plot_name="amino_acid_relative_counts",
                         count_dict=None,
                         scale=False,
                         lib_size_dict=None,
                         combine_sum=False, combine_weighted=False,
                         png_dir=png_dir, svg_dir=svg_dir):

        if count_dict is None:
            count_dict = self.amino_acid_df_dict

        if lib_size_dict is None:
            lib_size_dict = self.lib_size_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_heatmap_grid(title_prefix=title,
                               amino_acid_df_dict=count_dict,
                               scale=scale,
                               lib_size_dict=lib_size_dict,
                               combine_sum=combine_sum,
                               combine_weighted=combine_weighted,
                               png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_tabbed_heatmap_plot(self, group_count_dict,
                                plot_name="amino_acid_relative_counts",
                                scale=False, lib_size_dict_dict=None,
                                combine_sum=False, combine_weighted=False,
                                png_dir=png_dir, svg_dir=svg_dir):

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_tabbed_heatmap_grid(title_prefix=title,
                                      group_amino_acid_df_dict_dict=group_count_dict,
                                      scale=scale,
                                      lib_size_dict_dict=lib_size_dict_dict,
                                      combine_sum=combine_sum,
                                      combine_weighted=combine_weighted,
                                      png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_frame_barplot(self, plot_name="global_frame_preference",
                          count_dict=None, frame_stats_df_dict=None,
                          color_dict=None, lib_size_dict=None,
                          combine_sum=False, combine_weighted=False, combine_color=None,
                          png_dir=png_dir, svg_dir=svg_dir):

        if count_dict is None:
            count_dict = self.frame_count_start_dict

        if color_dict is None:
            color_dict = self.colors_dict

        if frame_stats_df_dict is None and not combine_sum and not combine_weighted:
            frame_stats_df_dict = self.frame_stats_df_dict

        if lib_size_dict is None:
            lib_size_dict = self.lib_size_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_frame_barplots(title_prefix=title,
                                 frame_df_dict=count_dict,
                                 frame_stats_df_dict=frame_stats_df_dict,
                                 lib_size_dict=lib_size_dict,
                                 color_dict=color_dict,
                                 combine_sum=combine_sum,
                                 combine_weighted=combine_weighted,
                                 combine_color=combine_color,
                                 png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_tabbed_frame_barplots(self, group_count_dict, group_frame_stats_df_dict,
                                  plot_name="global_frame_preference",
                                  color_dict=None, lib_size_dict_dict=None,
                                  combine_sum=False, combine_weighted=False, combine_color=None,
                                  png_dir=png_dir, svg_dir=svg_dir):

        if color_dict is None:
            color_dict = self.colors_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_tabbed_frame_barplots(title=title,
                                        group_frame_df_dict_dict=group_count_dict,
                                        group_frame_stats_df_dict_dict=group_frame_stats_df_dict,
                                        color_dict=color_dict,
                                        lib_size_dict_dict=lib_size_dict_dict,
                                        combine_sum=combine_sum,
                                        combine_weighted=combine_weighted,
                                        combine_color=combine_color,
                                        png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_fft_plot(self, plot_name="Periodicity-fast_Fourier_transform", region="",
                     count_dict=None,
                     color_dict=colors_dict,
                     lib_size_dict=None,
                     combine_sum=False, combine_weighted=False, combine_color=None,
                     png_dir=png_dir, svg_dir=svg_dir):

        if count_dict is None:
            if region == FivePSeqCounts.TERM:
                count_dict = self.fft_signal_term_dict
            else:
                count_dict = self.fft_signal_start_dict

        if lib_size_dict is None:
            lib_size_dict = self.lib_size_dict

        if color_dict is None:
            color_dict = self.colors_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name + "-" + region

        p = bokeh_fft_plot(title=title,
                           align_region=region,
                           signal_series_dict=count_dict,
                           color_dict=color_dict,
                           lib_size_dict=lib_size_dict,
                           combine_sum=combine_sum,
                           combine_weighted=combine_weighted,
                           combine_color=combine_color,
                           png_dir=png_dir, svg_dir=svg_dir)
        return p

    def get_tabbed_fft_plot(self, group_count_dict,
                            plot_name="Periodicity-fast_Fourier_transform",
                            region="",
                            color_dict=colors_dict,
                            scale=False, lib_size_dict_dict=None,
                            combine_sum=False, combine_weighted=False, combine_color=None,
                            png_dir=png_dir, svg_dir=svg_dir):

        if color_dict is None:
            color_dict = self.colors_dict

        if png_dir is None:
            png_dir = self.png_dir
        if svg_dir is None:
            svg_dir = self.svg_dir

        title = plot_name + "-" + region
        if scale:
            title += "_" + "scaled"
        else:
            title += "_" + "raw"

        if combine_weighted:
            title += "_" + "combined_weighted"
        elif combine_sum:
            title += "_" + "combined_summed"

        p = bokeh_tabbed_fft_plot(title=title,
                                  align_region=region,
                                  group_signal_series_dict_dict=group_count_dict,
                                  color_dict=color_dict,
                                  lib_size_dict_dict=lib_size_dict_dict,
                                  combine_sum=combine_sum,
                                  combine_weighted=combine_weighted,
                                  combine_color=combine_color,
                                  png_dir=png_dir, svg_dir=svg_dir)
        return p
