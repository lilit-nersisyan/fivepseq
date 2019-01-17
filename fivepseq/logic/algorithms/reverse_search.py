import os

import pandas as pd
import numpy as np
from bokeh.io import output_file, show
from bokeh.plotting import figure

from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from viz.bokeh_plots import bokeh_scatter_plot, bokeh_transcript_scatter_plot

dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"

transcript_assembly = pd.read_csv(
    os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "transcript_assembly.txt"),
    sep="\t")

transcript_count_full_Hela_rep1_dict = {
    "Hela-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_FULL_LENGTH.txt")),
    "HelaCHX-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_FULL_LENGTH.txt")),
    "HelaFrag-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_FULL_LENGTH.txt"))
}

HelaCHX_diff = [[]] * len(transcript_count_full_Hela_rep1_dict.get("HelaCHX-rep1"))
HelaUnt_diff = [[]] * len(transcript_count_full_Hela_rep1_dict.get("Hela-rep1"))

count_vector_list_frag = transcript_count_full_Hela_rep1_dict.get("HelaFrag-rep1")
count_vector_list_unt = transcript_count_full_Hela_rep1_dict.get("Hela-rep1")
count_vector_list_chx = transcript_count_full_Hela_rep1_dict.get("HelaCHX-rep1")

for i in range(0, len(transcript_assembly)):
    count_vector_unt = count_vector_list_unt[i]
    count_vector_chx = count_vector_list_chx[i]
    count_vector_frag = count_vector_list_frag[i]

    diff_vector_unt = np.array(count_vector_unt) - np.array(count_vector_frag)
    diff_vector_unt[diff_vector_unt < 0] = 0
    HelaUnt_diff[i] = list(diff_vector_unt)

    diff_vector_chx = np.array(count_vector_chx) - np.array(count_vector_frag)
    diff_vector_chx[diff_vector_chx < 0] = 0
    HelaCHX_diff[i] = list(diff_vector_chx)

chx_non_empty_ind = []
unt_non_empty_ind = []
for i in range(0, len(transcript_assembly)):
    if sum(HelaCHX_diff[i]) > 0 and sum(HelaUnt_diff[i]) > 0:
        chx_non_empty_ind.append(i)
        unt_non_empty_ind.append(i)

diff_count_dict = {
    "HelaCHX-rep1": HelaCHX_diff,
    "Hela-rep1": HelaUnt_diff
}

diff_color_dict = {
    "HelaCHX-rep1": "green",
    "Hela-rep1": "blue"
}

for i in range(0, len(transcript_assembly)):
    np.fft(HelaCHX_diff[i])

title = "HelaCHX-rep1_fft"
output_file(title + ".html")
x = np.fft.fft(HelaCHX_diff[6])
p = figure(title=title, x_axis_label="frequency", y_axis_label="value")
line = p.line(list(range(100)), abs(x[range(100)]), line_color="red")
show(p)

bokeh_transcript_scatter_plot("Hela-rep1_diffs", diff_count_dict, transcript_assembly, diff_color_dict,
                              FivePSeqCounts.TERM, 500, index_filter=chx_non_empty_ind[slice(0, 100)], min_count=50)
