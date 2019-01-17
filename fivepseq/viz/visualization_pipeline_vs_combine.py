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
            frame_count_START_combined.loc[:,('F0', 'F1', 'F2')] += frame_start.loc[:,('F0', 'F1', 'F2')]
            frame_count_TERM_combined.loc[:,('F0', 'F1', 'F2')] += frame_term.loc[:,('F0', 'F1', 'F2')]
            amino_acid_df_combined += amino_acid_df

    
    # combined plots

    colors_dict = {group + "_combined" : cl.to_numeric(cl.scales['9']['qual']['Set3'][3])}

    p_scatter_term = bokeh_scatter_plot(group + "_combined_" + org + "_term", FivePSeqCounts.TERM,
                                        {group + "_combined" : meta_count_term_combined},
                                        colors_dict)
    p_scatter_start = bokeh_scatter_plot(group + "_combined_" + org + "_start", FivePSeqCounts.START,
                                         {group + "_combined": meta_count_start_combined},
                                         colors_dict)
    p_triangle = bokeh_triangle_plot(group + "_combined_" + org + "_triangle",
                                     {group + "_combined": frame_count_TERM_combined}, colors_dict)
    p_triangle_START = bokeh_triangle_plot(group + "_combined_" + org + "_triangle_START",
                                           {group + "_combined": frame_count_START_combined},
                                           colors_dict)
    p_aa_heatmaps = bokeh_heatmap_grid(group + "_combined_" + org + "_amino_acid_pauses",
                                       {group + "_combined": amino_acid_df_combined})

    title = group + "_" + org
    bokeh_composite(group + "_combined_" + org, [p_scatter_start, p_scatter_term,
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
