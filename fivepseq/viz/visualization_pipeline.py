import os

import colorlover as cl
import pandas as pd
from fivepseq.util.writers import FivePSeqOut

from fivepseq.logic.structures.fivepseq_counts import CountManager
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts
from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_composite, bokeh_heatmap_grid, \
    bokeh_frame_barplots

FILTER_TOP_POPULATED = "Filter top populated transcripts"
FILTER_CANONICAL_TRANSCRIPTS = "Filter canonical transcripts"


def plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots, subgroup=""):
    print "organism: " + org
    transcript_count_full_dict = {}
    meta_count_term_dict = {}
    meta_count_start_dict = {}
    amino_acid_df_dict = {}
    frame_count_dict = {}
    frame_count_dict_START = {}
    for sample in samples:
        transcript_count_full_dict.update({sample: CountManager.read_counts_as_list(
            os.path.join(dir_5pseq_group, sample + "-" + org, "counts_FULL_LENGTH.txt"))})

        meta_count_term_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "meta_counts_TERM.txt"),
            sep="\t", header=None, names=["D", "C"])})

        meta_count_start_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "meta_counts_START.txt"),
            sep="\t", header=None, names=["D", "C"])})

        frame_count_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "frame_counts_TERM.txt"),
            sep="\t")})

        frame_count_dict_START.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "frame_counts_START.txt"),
            sep="\t")})

        if os.path.exists(os.path.join(dir_5pseq_group, sample + "-" + org, FivePSeqOut.AMINO_ACID_PAUSES_FILE)):
            amino_acid_df_dict.update({sample: pd.read_csv(
                os.path.join(dir_5pseq_group, sample + "-" + org, "amino_acid_pauses.txt"),
                sep="\t", header=0, index_col=0
            )})
        else:
            amino_acid_df_dict.update({sample: None})

    large_colors_list = (cl.to_numeric(cl.scales['8']['qual']['Paired']) +
                         cl.to_numeric(cl.scales['8']['qual']['Set1']) +
                         cl.to_numeric(cl.scales['5']['qual']['Set3']))

    colors_dict = dict(
        zip(transcript_count_full_dict.keys(),
            large_colors_list[0:len(samples)]))

    p_scatter_term = bokeh_scatter_plot(group + "_" + org + "_term", FivePSeqCounts.TERM, meta_count_term_dict,
                                        colors_dict)
    p_scatter_start = bokeh_scatter_plot(group + "_" + org + "_start", FivePSeqCounts.START, meta_count_start_dict,
                                         colors_dict)
    p_triangle = bokeh_triangle_plot(group + "_" + org + "_triangle", frame_count_dict, colors_dict)
    p_triangle_START = bokeh_triangle_plot(group + "_" + org + "_triangle_START", frame_count_dict_START, colors_dict)
    # print "Plotting meta scatter plot ..."
    # p_normalized_meta_scatter_plot = bokeh_normalized_meta_scatter_plot(group + "_" + org + "_norm_meta_counts",
    #                                                                    transcript_count_full_dict, colors_dict,
    #                                                                    x_max=1000, span_size=100)
    p_aa_heatmaps = bokeh_heatmap_grid(group + "_" + org + "_amino_acid_pauses", amino_acid_df_dict)
    p_frame_barplots = bokeh_frame_barplots(group + "_" + org + "_frame_histograms", frame_count_dict, colors_dict)
    p_frame_barplots_start = bokeh_frame_barplots(group + "_" + org + "_frame_START_histograms",
                                                  frame_count_dict_START, colors_dict)

    ######          Combined plots              ###############

    meta_count_start_combined = pd.DataFrame()
    meta_count_term_combined = pd.DataFrame()
    frame_count_START_combined = pd.DataFrame()
    frame_count_TERM_combined = pd.DataFrame()
    amino_acid_df_combined = pd.DataFrame()
    for key in meta_count_start_dict.keys():
        df_start = meta_count_start_dict.get(key)
        df_term = meta_count_term_dict.get(key)
        frame_start = frame_count_dict_START.get(key)
        frame_term = frame_count_dict.get(key)
        amino_acid_df = amino_acid_df_dict.get(key)
        if len(meta_count_start_combined) == 0:
            meta_count_start_combined = df_start.copy()
            meta_count_term_combined = df_term.copy()
            frame_count_START_combined = frame_start.copy()
            frame_count_TERM_combined = frame_term.copy()
            amino_acid_df_combined = amino_acid_df.copy()
        else:
            meta_count_start_combined.C += df_start.C
            meta_count_term_combined.C += df_term.C
            frame_count_START_combined.loc[:, ('F0', 'F1', 'F2')] += frame_start.loc[:, ('F0', 'F1', 'F2')]
            frame_count_TERM_combined.loc[:, ('F0', 'F1', 'F2')] += frame_term.loc[:, ('F0', 'F1', 'F2')]
            amino_acid_df_combined += amino_acid_df

    colors_dict = {group + "_combined": cl.to_numeric(cl.scales['9']['qual']['Set3'])[3]}

    p_scatter_term_combined = bokeh_scatter_plot(group + "_combined_" + org + "_term", FivePSeqCounts.TERM,
                                                 {group + "_combined": meta_count_term_combined},
                                                 colors_dict)
    p_scatter_start_combined = bokeh_scatter_plot(group + "_combined_" + org + "_start", FivePSeqCounts.START,
                                                  {group + "_combined": meta_count_start_combined},
                                                  colors_dict)
    p_triangle_combined = bokeh_triangle_plot(group + "_combined_" + org + "_triangle",
                                              {group + "_combined": frame_count_TERM_combined}, colors_dict)
    p_triangle_START_combined = bokeh_triangle_plot(group + "_combined_" + org + "_triangle_START",
                                                    {group + "_combined": frame_count_START_combined},
                                                    colors_dict)
    p_aa_heatmaps_combined = bokeh_heatmap_grid(group + "_combined_" + org + "_amino_acid_pauses",
                                                {group + "_combined": amino_acid_df_combined})

    p_frame_barplots_combined = bokeh_frame_barplots(group + "_combined_" + org + "_frame_histograms",
                                                     {group + "_combined": frame_count_TERM_combined},
                                                     colors_dict)
    p_frame_barplots_start_combined = bokeh_frame_barplots(group + "_combined_" + org + "_frame_START_histograms",
                                                           {group + "_combined": frame_count_START_combined},
                                                           colors_dict)
    if subgroup == "":
        title = group + "_" + org
    else:
        title = group + "-" + subgroup + "_" + org
    bokeh_composite(group + "_" + org, [p_scatter_start, p_scatter_term,
                                        p_triangle, p_triangle_START,
                                        p_aa_heatmaps, None, None, None,
                                        p_frame_barplots, None, None, None,
                                        p_frame_barplots_start, None, None, None,
                                        p_scatter_start_combined, p_scatter_term_combined,
                                        p_triangle_combined, p_triangle_START_combined,
                                        p_aa_heatmaps_combined, None, None, None,
                                        p_frame_barplots_combined, None, None, None,
                                        p_frame_barplots_start_combined, None, None, None],
                    os.path.join(dir_5pseq_plots, title + ".html"), 4)


def plot_single_sample(sample_name, title, count_dir, plot_dir, span_size=500,
                       transcript_index=None, filter=None):
    fivepseqout = FivePSeqOut(count_dir)
    meta_count_term = CountManager.read_meta_counts(
        fivepseqout.get_file_path(FivePSeqOut.META_COUNT_TERM_FILE))
    meta_count_start = CountManager.read_meta_counts(
        fivepseqout.get_file_path(FivePSeqOut.META_COUNT_START_FILE))

    try:
        amino_acid_df = CountManager.read_amino_acid_df(
            fivepseqout.get_file_path(FivePSeqOut.AMINO_ACID_PAUSES_FILE))
    except:
        amino_acid_df = None
    frame_count_term = CountManager.read_frame_counts(
        fivepseqout.get_file_path(FivePSeqOut.FRAME_COUNTS_TERM_FILE))
    frame_count_start = CountManager.read_frame_counts(
        fivepseqout.get_file_path(FivePSeqOut.FRAME_COUNTS_START_FILE))

    # count_vector_list = CountManager.read_counts_as_list(
    #    fivepseqout.get_file_path(FivePSeqOut.COUNT_FULL_FILE))
    count_vector_list_term = CountManager.read_counts_as_list(
        fivepseqout.get_file_path(FivePSeqOut.COUNT_TERM_FILE))
    count_vector_list_start = CountManager.read_counts_as_list(
        fivepseqout.get_file_path(FivePSeqOut.COUNT_START_FILE))

    if filter is not None:
        if filter == FILTER_TOP_POPULATED:
            transcript_index = CountManager.top_populated_count_vector_indices(
                count_vector_list_term, span_size, 1000)
        elif filter == FILTER_CANONICAL_TRANSCRIPTS:
            transcript_index = CountManager.canonical_transcript_indices(count_dir)

    if transcript_index is not None:
        frame_count_term = frame_count_term.iloc[transcript_index,]
        frame_count_start = frame_count_start.iloc[transcript_index,]
        count_vector_list_term = [count_vector_list_term[i] for i in transcript_index]
        count_vector_list_start = [count_vector_list_start[i] for i in transcript_index]
        meta_count_term = CountManager.count_vector_to_df(
            CountManager.compute_meta_counts(count_vector_list_term),
            FivePSeqCounts.TERM, span_size)
        meta_count_start = CountManager.count_vector_to_df(
            CountManager.compute_meta_counts(count_vector_list_start),
            FivePSeqCounts.TERM, span_size)

        # TODO amino acids pauses not subsettable

    color = cl.to_numeric(cl.scales['5']['qual']['Set3'])[0]

    p_scatter_term = bokeh_scatter_plot(sample_name + "_term", FivePSeqCounts.TERM,
                                        {sample_name: meta_count_term},
                                        {sample_name: color})
    p_scatter_start = bokeh_scatter_plot(sample_name + "_start", FivePSeqCounts.START,
                                         {sample_name: meta_count_start},
                                         {sample_name: color})
    p_triangle = bokeh_triangle_plot(sample_name + "_triangle",

                                     {sample_name: frame_count_term},
                                     {sample_name: color})

    p_triangle_START = bokeh_triangle_plot(sample_name + "_triangle_START",
                                           {sample_name: frame_count_start},
                                           {sample_name: color})

    p_aa_heatmap = bokeh_heatmap_grid(sample_name + "_amino_acid_pauses",
                                      {sample_name: amino_acid_df})
    p_frame_barplots_term = bokeh_frame_barplots(sample_name + "_frame_histograms",
                                                 {sample_name: frame_count_term},
                                                 {sample_name: color})
    p_frame_barplots_start = bokeh_frame_barplots(sample_name + "_frame_START_histograms",
                                                  {sample_name: frame_count_start},
                                                  {sample_name: color})

    bokeh_composite(title, [p_scatter_start, p_scatter_term,
                            p_triangle, p_triangle_START,
                            p_aa_heatmap, None,
                            p_frame_barplots_term, p_frame_barplots_start],
                    os.path.join(plot_dir, title + ".html"), 2)


def plot_hs_cell_lines(cell_line, samples, dir_5pseq_group, dir_5pseq_plots, subgroup=""):
    print "cell line: " + cell_line
    transcript_count_full_dict = {}
    meta_count_term_dict = {}
    meta_count_start_dict = {}
    amino_acid_df_dict = {}
    frame_count_dict = {}
    frame_count_dict_START = {}
    for sample in samples:
        transcript_count_full_dict.update({sample: CountManager.read_counts_as_list(
            os.path.join(dir_5pseq_group + sample, "counts_FULL_LENGTH.txt"))})

        meta_count_term_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group + sample, "meta_counts_TERM.txt"),
            sep="\t", header=None, names=["D", "C"])})

        meta_count_start_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group + sample, "meta_counts_START.txt"),
            sep="\t", header=None, names=["D", "C"])})

        frame_count_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group + sample, "frame_counts_TERM.txt"),
            sep="\t")})

        frame_count_dict_START.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group + sample, "frame_counts_START.txt"),
            sep="\t")})

        amino_acid_df_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group + sample, "amino_acid_pauses.txt"),
            sep="\t", header=0, index_col=0
        )})

    large_colors_list = (cl.to_numeric(cl.scales['8']['qual']['Paired']) +
                         cl.to_numeric(cl.scales['8']['qual']['Set1']) +
                         cl.to_numeric(cl.scales['5']['qual']['Set3']))

    colors_dict = dict(
        zip(transcript_count_full_dict.keys(),
            large_colors_list[0:len(samples)]))

    p_scatter_term = bokeh_scatter_plot(cell_line + "_term", FivePSeqCounts.TERM, meta_count_term_dict,
                                        colors_dict)
    p_scatter_start = bokeh_scatter_plot(cell_line + "_start", FivePSeqCounts.START, meta_count_start_dict,
                                         colors_dict)
    p_triangle = bokeh_triangle_plot(cell_line + "_triangle", frame_count_dict, colors_dict)
    p_triangle_START = bokeh_triangle_plot(cell_line + "_triangle_START", frame_count_dict_START, colors_dict)
    # print "Plotting meta scatter plot ..."
    # p_normalized_meta_scatter_plot = bokeh_normalized_meta_scatter_plot(group + "_" + org + "_norm_meta_counts",
    #                                                                    transcript_count_full_dict, colors_dict,
    #                                                                    x_max=1000, span_size=100)
    p_aa_heatmaps = bokeh_heatmap_grid(cell_line + "_amino_acid_pauses", amino_acid_df_dict)
    p_frame_barplots = bokeh_frame_barplots(cell_line + "_frame_histograms", frame_count_dict, colors_dict)
    p_frame_barplots_start = bokeh_frame_barplots(cell_line + "_frame_START_histograms",
                                                  frame_count_dict_START, colors_dict)

    ######          Combined plots              ###############

    meta_count_start_combined = pd.DataFrame()
    meta_count_term_combined = pd.DataFrame()
    frame_count_START_combined = pd.DataFrame()
    frame_count_TERM_combined = pd.DataFrame()
    amino_acid_df_combined = pd.DataFrame()
    for key in meta_count_start_dict.keys():
        if "Frag" not in key:
            df_start = meta_count_start_dict.get(key)
            df_term = meta_count_term_dict.get(key)
            frame_start = frame_count_dict_START.get(key)
            frame_term = frame_count_dict.get(key)
            amino_acid_df = amino_acid_df_dict.get(key)
            if len(meta_count_start_combined) == 0:
                meta_count_start_combined = df_start.copy()
                meta_count_term_combined = df_term.copy()
                frame_count_START_combined = frame_start.copy()
                frame_count_TERM_combined = frame_term.copy()
                amino_acid_df_combined = amino_acid_df.copy()
            else:
                meta_count_start_combined.C += df_start.C
                meta_count_term_combined.C += df_term.C
                frame_count_START_combined.loc[:, ('F0', 'F1', 'F2')] += frame_start.loc[:, ('F0', 'F1', 'F2')]
                frame_count_TERM_combined.loc[:, ('F0', 'F1', 'F2')] += frame_term.loc[:, ('F0', 'F1', 'F2')]
                amino_acid_df_combined += amino_acid_df

    colors_dict = {cell_line + "_combined": cl.to_numeric(cl.scales['9']['qual']['Set3'])[3]}

    p_scatter_term_combined = bokeh_scatter_plot(cell_line + "_combined_" + "_term", FivePSeqCounts.TERM,
                                                 {cell_line + "_combined": meta_count_term_combined},
                                                 colors_dict)
    p_scatter_start_combined = bokeh_scatter_plot(cell_line + "_combined_" + "_start", FivePSeqCounts.START,
                                                  {cell_line + "_combined": meta_count_start_combined},
                                                  colors_dict)
    p_triangle_combined = bokeh_triangle_plot(cell_line + "_combined_" + "_triangle",
                                              {cell_line + "_combined": frame_count_TERM_combined}, colors_dict)
    p_triangle_START_combined = bokeh_triangle_plot(cell_line + "_combined_" + "_triangle_START",
                                                    {cell_line + "_combined": frame_count_START_combined},
                                                    colors_dict)
    p_aa_heatmaps_combined = bokeh_heatmap_grid(cell_line + "_combined_" + "_amino_acid_pauses",
                                                {cell_line + "_combined": amino_acid_df_combined})

    p_frame_barplots_combined = bokeh_frame_barplots(cell_line + "_combined_" + "_frame_histograms",
                                                     {cell_line + "_combined": frame_count_TERM_combined},
                                                     colors_dict)
    p_frame_barplots_start_combined = bokeh_frame_barplots(cell_line + "_combined_" + "_frame_START_histograms",
                                                           {cell_line + "_combined": frame_count_START_combined},
                                                           colors_dict)
    if subgroup == "":
        title = cell_line
    else:
        title = cell_line + "-" + subgroup
    bokeh_composite(cell_line, [p_scatter_start, p_scatter_term,
                                p_triangle, p_triangle_START,
                                p_aa_heatmaps, None, None, None,
                                p_frame_barplots, None, None, None,
                                p_frame_barplots_start, None, None, None,
                                p_scatter_start_combined, p_scatter_term_combined,
                                p_triangle_combined, p_triangle_START_combined,
                                p_aa_heatmaps_combined, None, None, None,
                                p_frame_barplots_combined, None, None, None,
                                p_frame_barplots_start_combined, None, None, None],
                    os.path.join(dir_5pseq_plots, title + ".html"), 4)
