import os
from sets import Set

import pandas as pd
import numpy as np
import colorlover as cl
from fivepseq.logic.structures.counts import FivePSeqCounts

from fivepseq.viz.bokeh_plots import bokeh_scatter_plot, bokeh_triangle_plot, bokeh_transcript_scatter_plot, \
    bokeh_composite, bokeh_normalized_meta_scatter_plot

from fivepseq.logic.structures.fivepseq_counts import CountManager

dir_5pseq_microbiome = "/proj/sllstore2017018/lilit/5pseq_microbiome"
meta_count_dict_lpla = {
    "lpla_untreated": pd.read_csv(
        os.path.join(dir_5pseq_microbiome, "fivepseq_lpla_untreated", "meta_counts_TERM.txt"),
        sep="\t"),
    "lpla_fragmented": pd.read_csv(
        os.path.join(dir_5pseq_microbiome, "fivepseq_l_pla_fragmented", "meta_counts_TERM.txt"),
        sep="\t")
}
colors_dict_lpla = dict(zip(meta_count_dict_lpla.keys(),
                            np.array(cl.scales['3']['div']['RdGy'])[[0, 2]]))

dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"

meta_count_Hela_rep1_term_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_TERM.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_TERM.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_TERM.txt"),
        sep="\t", header=None, names=["D", "C"]),
}

meta_count_Hela_rep1_start_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_START.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_START.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_START.txt"),
        sep="\t", header=None, names=["D", "C"]),
}

colors_dict_Hela = dict(zip(meta_count_Hela_rep1_term_dict.keys(), cl.to_numeric(cl.scales['6']['qual']['Set1'])))

frame_count_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "frame_counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "frame_counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "frame_counts_TERM.txt"),
        sep="\t"),
}

transcript_count_term_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_TERM.txt"),
        sep="\t"),
}

transcript_count_start_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_START.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_START.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_START.txt"),
        sep="\t"),
}

transcript_count_full_Hela_rep1_dict = {
    "Hela-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_FULL_LENGTH.txt")),
    "HelaCHX-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_FULL_LENGTH.txt")),
    "HelaFrag-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_FULL_LENGTH.txt"))
}

#########################################
#          filter transcripts           #
#########################################

transcript_assembly = pd.read_csv(
    os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "transcript_assembly.txt"),
    sep="\t")

filtered_index = []
gene_set = Set()
unique_transcripts = []
for i in range(0, len(transcript_assembly)):
    gene = transcript_assembly.gene[i].split("gene:")[1]
    gene_set.add(gene)

gene_transcript_dict = {}
for i in range(0, len(transcript_assembly)):
    gene = transcript_assembly.gene[i].split("gene:")[1]
    transcript = transcript_assembly.ID[i].split("transcript:")[1]
    if gene_transcript_dict.has_key(gene):
        gene_transcript_dict.get(gene).append(transcript)
    else:
        gene_transcript_dict.update({gene: [transcript]})

for i in range(0, len(transcript_assembly)):
    gene = transcript_assembly.gene[i].split("gene:")[1]
    if len(gene_transcript_dict.get(gene)) == 1:
        unique_transcripts.append(i)

p_scatter_term = bokeh_scatter_plot("Hela-rep1_term", FivePSeqCounts.TERM, meta_count_Hela_rep1_term_dict,
                                    colors_dict_Hela)
p_scatter_start = bokeh_scatter_plot("Hela-rep1_start", FivePSeqCounts.START, meta_count_Hela_rep1_start_dict,
                                     colors_dict_Hela)
p_triangle = bokeh_triangle_plot("Hela-rep1_triangle", frame_count_Hela_rep1_dict, colors_dict_Hela,
                                 transcript_index=unique_transcripts)
p_transcripts = bokeh_transcript_scatter_plot("Hela-rep1_transcript_counts",
                                              transcript_count_full_Hela_rep1_dict,
                                              transcript_assembly,
                                              colors_dict_Hela,
                                              FivePSeqCounts.TERM, 500, unique_transcripts, 300, 500)
bokeh_composite("Hela-rep1", [p_scatter_start, p_scatter_term, p_triangle, p_transcripts], 3)

#########################################
#          zymo
#########################################

dir_5pseq_zymo = "/proj/sllstore2017018/lilit/vagaga2/zymo/zymo_fivepseq_best"

def plot_samples_for_organism(org):
    msc1_dir = "MCS1-" + org + "_fivepseq"
    msc2_dir = "MCS2-" + org + "_fivepseq"
    transcript_count_full_zymo_b_sub_dict = {
        "MCS1": CountManager.read_counts_as_list(
            os.path.join(dir_5pseq_zymo, msc1_dir, "counts_FULL_LENGTH.txt")),
        "MCS2": CountManager.read_counts_as_list(
            os.path.join(dir_5pseq_zymo, msc2_dir, "counts_FULL_LENGTH.txt"))
    }

    meta_count_zymo_b_sub_term_dict = {
        "MCS1": pd.read_csv(
            os.path.join(dir_5pseq_zymo, msc1_dir, "meta_counts_TERM.txt"),
            sep="\t", header=None, names=["D", "C"]),
        "MCS2": pd.read_csv(
            os.path.join(dir_5pseq_zymo, msc2_dir, "meta_counts_TERM.txt"),
            sep="\t", header=None, names=["D", "C"])
    }

    meta_count_zymo_b_sub_start_dict = {
        "MCS1": pd.read_csv(
            os.path.join(dir_5pseq_zymo, msc1_dir, "meta_counts_START.txt"),
            sep="\t", header=None, names=["D", "C"]),
        "MCS2": pd.read_csv(
            os.path.join(dir_5pseq_zymo, msc2_dir, "meta_counts_START.txt"),
            sep="\t", header=None, names=["D", "C"])
    }

    frame_count_zymo_b_sub_dict = {
        "MCS1": pd.read_csv(
            os.path.join(dir_5pseq_zymo, msc1_dir, "frame_counts_TERM.txt"),
            sep="\t"),
        "MCS2": pd.read_csv(
            os.path.join(dir_5pseq_zymo, msc2_dir, "frame_counts_TERM.txt"),
            sep="\t"),
    }

    colors_dict_zymo = dict(
        zip(transcript_count_full_zymo_b_sub_dict.keys(), cl.to_numeric(cl.scales['6']['qual']['Set1'])))

    p_scatter_term = bokeh_scatter_plot("zymo_" + org + "_term", FivePSeqCounts.TERM, meta_count_zymo_b_sub_term_dict,
                                        colors_dict_zymo)
    p_scatter_start = bokeh_scatter_plot("zymo_" + org + "_start", FivePSeqCounts.START, meta_count_zymo_b_sub_start_dict,
                                         colors_dict_zymo)
    p_triangle = bokeh_triangle_plot("zymo_" + org + "_triangle", frame_count_zymo_b_sub_dict, colors_dict_zymo)
    p_normalized_meta_scatter_plot = bokeh_normalized_meta_scatter_plot("zymo_" + org + "_norm_meta_counts", transcript_count_full_zymo_b_sub_dict, colors_dict_zymo, 200)
    bokeh_composite("zymo_" + org, [p_scatter_start, p_scatter_term, p_triangle, p_normalized_meta_scatter_plot], 2)

organisms = (
    "Bacillus_subtilis",
    "Cryptococcus_neoformans",
    "Enterococcus_faecalis",
    "Escherichia_coli",
    "Lactobacillus_fermentum",
    "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa",
    "Saccharomyces_cerevisiae",
    "Salmonella_enterica",
    "Staphylococcus_aureus"
)
for org in organisms:
    print org
    plot_samples_for_organism(org)
