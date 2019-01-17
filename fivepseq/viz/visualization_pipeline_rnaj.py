import os
from sets import Set

import pandas as pd
import numpy as np
import colorlover as cl
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts

from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_transcript_scatter_plot, \
    bokeh_composite, bokeh_normalized_meta_scatter_plot, bokeh_heatmap_grid
from bokeh.palettes import Category20b
from fivepseq.logic.structures.fivepseq_counts import CountManager

#########################################
#          zymo
#########################################
group = "rnaj"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

def plot_samples_for_organism(organism, samples, group, subgroup):
    print "organism: " + organism
    transcript_count_full_dict = {}
    meta_count_term_dict = {}
    meta_count_start_dict = {}
    frame_count_dict = {}
    frame_count_dict_START = {}
    amino_acid_df_dict = {}
    for sample in samples:
        transcript_count_full_dict.update({sample: CountManager.read_counts_as_list(
            os.path.join(dir_5pseq_group, sample + "-" + organism, "counts_FULL_LENGTH.txt"))})

        meta_count_term_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + organism, "meta_counts_TERM.txt"),
            sep="\t", header=None, names=["D", "C"])})

        meta_count_start_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + organism, "meta_counts_START.txt"),
            sep="\t", header=None, names=["D", "C"])})

        frame_count_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + organism, "frame_counts_TERM.txt"),
            sep="\t")})

        frame_count_dict_START.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "frame_counts_START.txt"),
            sep="\t")})

        amino_acid_df_dict.update({sample: pd.read_csv(
            os.path.join(dir_5pseq_group, sample + "-" + org, "amino_acid_pauses.txt"),
            sep="\t", header=0, index_col=0
        )})

    mycolors = list(Category20b[20][i] for i in
                    list(range(0, 3)) + list(range(4, 7)) + list(range(8, 11)) + list(range(12, 15)) + list(
                        range(16, 19)))
    colors_dict = dict(
        zip(transcript_count_full_dict.keys(),
            cl.to_numeric(
                cl.scales['3']['seq']['Purples'] + cl.scales['3']['seq']['Greens'] +
                cl.scales['3']['seq']['Reds'] + cl.scales['3']['seq']['Blues'] +
                cl.scales['3']['seq']['Greys'])))

    print "Plotting terminal meta counts ..."
    p_scatter_term = bokeh_scatter_plot(group + "_" + organism + "_term", FivePSeqCounts.TERM, meta_count_term_dict,
                                        colors_dict)
    print "Plotting start meta counts ..."
    p_scatter_start = bokeh_scatter_plot(group + "_" + organism + "_start", FivePSeqCounts.START, meta_count_start_dict,
                                         colors_dict)
    print "Plotting triangle plot ..."
    p_triangle = bokeh_triangle_plot(group + "_" + organism + "_triangle", frame_count_dict, colors_dict)
    p_triangle_START = bokeh_triangle_plot(group + "_" + org + "_triangle_START", frame_count_dict_START, colors_dict)
    #print "Plotting meta scatter plot ..."
    #p_normalized_meta_scatter_plot = bokeh_normalized_meta_scatter_plot(group + "_" + org + "_norm_meta_counts",
    #                                                                    transcript_count_full_dict, colors_dict,
    #                                                                    x_max=1000, span_size=100)
    p_aa_heatmaps = bokeh_heatmap_grid(group + "_" + org + "_amino_acid_pauses", amino_acid_df_dict)

    title = group + "-" + subgroup + "_" + org
    bokeh_composite(group + "_" + org, [p_scatter_start, p_scatter_term,
                                        p_triangle, p_triangle_START, p_aa_heatmaps],
                    os.path.join(dir_5pseq_plots, title + ".html"), 4)



samples = (
    "rnja1delta_168trpC2-CAM1",
    "rnja1delta_168trpC2-CAM2",
    "rnja1delta_168trpC2-CAM3",
    "rnja1delta_168trpC2-Ctr1",
    "rnja1delta_168trpC2-Ctr2",
    "rnja1delta_168trpC2-Ctr3",
)


samples_wt = (
    "WT_168trpC2-CAM1",
    "WT_168trpC2-CAM2",
    "WT_168trpC2-CAM3",
    "WT_168trpC2-Ctr1",
    "WT_168trpC2-Ctr2",
    "WT_168trpC2-Ctr3",
    "WT_168trpC2-fragCtr1",
    "WT_168trpC2-fragCtr2",
    "WT_168trpC2-fragCtr3",
)





organisms = (
    "Bacillus_subtilis_subsp_subtilis_str_168",
)
for org in organisms:
    plot_samples_for_organism(org, samples, group, "rnaj")
    plot_samples_for_organism(org, samples_wt, group, "WT")
