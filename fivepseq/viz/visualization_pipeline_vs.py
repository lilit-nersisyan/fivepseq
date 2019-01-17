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
group = "vs"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"

if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "VS_106_S25",
    "VS_153_S26",
    "VS_156_S27",
    "VS_162_S28",
    "VS_189_S29",
    "VS_35_S22",
    "VS_64_S23",
    "VS_88_S24",
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

    colors_dict = dict(
        zip(transcript_count_full_dict.keys(),
            cl.to_numeric(cl.scales['8']['qual']['Set1'])))

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

    title = group + "_" + org
    bokeh_composite(group + "_" + org, [p_scatter_start, p_scatter_term,
                                        p_triangle, p_triangle_START, p_aa_heatmaps],
                    os.path.join(dir_5pseq_plots, title + ".html"), 4)


organisms = (
    "Alteromonas_australica",
    "Bacillus_subtilis",
    "Escherichia_coli",
    "Gardnerella_vaginalis",
    "Lactobacillus_amylophilus",
    "Lactobacillus_crispatus",
    "Lactobacillus_jensenii",
    "Mageeibacillus_indolicus",
    "Prevotella_denticola",
    "Prevotella_fusca",
    "Prevotella_sp_oral_taxon_299",
    "Ruminococcus_champanellensis",
    "Salmonella_enterica",
    "Serratia_marcescens",
    "Ureaplasma_parvum",
    "Xanthomonas_euvesicatoria",

)
for org in organisms:
    plot_samples_for_organism(org)
