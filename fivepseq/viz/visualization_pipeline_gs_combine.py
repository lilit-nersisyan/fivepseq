import os
from sets import Set

import pandas as pd
import numpy as np
import colorlover as cl
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts

from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_transcript_scatter_plot, \
    bokeh_composite, bokeh_normalized_meta_scatter_plot, bokeh_heatmap_grid

from fivepseq.logic.structures.fivepseq_counts import CountManager

#########################################
#          zymo
#########################################
group = "gs"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"

if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "2GS_X2_Antrum_TOT_S21",
    "2GS_X4_Fundus_P16_S11",
    "GS_X3_Antrum_P16_S8",
    "GS_X3_MID_P16_S5",
    "2GS_X2_Fundus_TOT_S19",
    "2GS_X4_Fundus_P8_S10",
    "GS_X3_Antrum_P8_S7",
    "GS_X3_MID_P8_S4",
    "2GS_X2_MID_TOT_S20",
    "2GS_X4_Fundus_S16_S12",
    "GS_X3_Antrum_S16_S9",
    "GS_X3_MID_S16_S6",
    "2GS_X4_Antrum_P16_S17",
    "2GS_X4_MID_P16_S14",
    "GS_X3_Fundus_P16_S2",
    "2GS_X4_Antrum_P8_S16",
    "2GS_X4_MID_P8_S13",
    "GS_X3_Fundus_P8_S1",
    "2GS_X4_Antrum_S16_S18",
    "2GS_X4_MID_S16_S15",
    "GS_X3_Fundus_S16_S3"
)


def plot_samples_for_organism(org):
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

        amino_acid_df_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "amino_acid_pauses.txt"),
            sep="\t", header=0, index_col=0
        )})

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

    colors_dict = dict(
        zip(transcript_count_full_dict.keys(),
            cl.to_numeric(cl.scales['8']['qual']['Paired']) +
            cl.to_numeric(cl.scales['8']['qual']['Set1']) +
            cl.to_numeric(cl.scales['5']['qual']['Set3'])))

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

    # combined plots

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

    title = group + "_" + org
    bokeh_composite(group + "_" + org, [p_scatter_start, p_scatter_term,
                                        p_triangle, p_triangle_START,
                                        p_aa_heatmaps, None, None, None,
                                        p_scatter_start_combined, p_scatter_term_combined,
                                        p_triangle_combined, p_triangle_START_combined,
                                        p_aa_heatmaps_combined],
                    os.path.join(dir_5pseq_plots, title + ".html"), 4)


organisms = (
    "Bacillus_anthracis",
    "Bacillus_cereus",
    "Bacillus_subtilis",
    "Bacillus_thuringiensis",
    "Bacillus_velezensis",
    "Escherichia_coli",
    "Klebsiella_pneumoniae",
    "Salmonella_enterica",
    "Serratia_marcescens"

)
for org in organisms:
    plot_samples_for_organism(org)
