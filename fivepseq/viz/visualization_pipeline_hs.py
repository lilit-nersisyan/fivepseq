import os

from fivepseq.logic.structures.fivepseq_counts import CountManager

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism, plot_hs_cell_lines, plot_single_sample, \
    FILTER_TOP_POPULATED, FILTER_CANONICAL_TRANSCRIPTS

dir_5pseq_plots = "/proj/sllstore2017018/lilit/5pseq_human/fivepseq_plots"

if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "-rep1",
    "-rep2",
    "CHX-rep1",
    "CHX-rep2",
    "Frag-rep1",
    "Frag-rep2",
)

cell_lines = {
    "Hela",
    "HEK293"
}

"""
for cell_line in cell_lines:
    dir_5pseq_group = "/proj/sllstore2017018/lilit/5pseq_human/fivepseq_" + cell_line
    dir_5pseq_plots = "/proj/sllstore2017018/lilit/5pseq_human/fivepseq_plots"
    #plot_hs_cell_lines(cell_line, samples, dir_5pseq_group, dir_5pseq_plots)
    for sample in samples:
        plot_single_sample(cell_line + sample, cell_line + sample + "_top_populated",
                           dir_5pseq_group + sample,
                           dir_5pseq_plots,
                           span_size=500,
                           filter=FILTER_TOP_POPULATED)

count_dir = "/proj/sllstore2017018/lilit/vagaga2/lreu/lreu_best_fivepseq/lreu_10000-Lactobacillus_reuteri"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/5pseq_human/fivepseq_plots"
sample = "lreu_10000-Lactobacillus_reuteri"
plot_single_sample(sample, sample + "_top_populated",
                   count_dir, dir_5pseq_plots, 100, filter=FILTER_TOP_POPULATED)
"""

sample = "fivepseq_HelaCHX-rep1_test"
plot_single_sample(sample, sample + "_canonical_transcripts",
                   "/proj/sllstore2017018/lilit/5pseq_human/fivepseq_HelaCHX-rep1_test",
                   dir_5pseq_plots, 500, filter=FILTER_CANONICAL_TRANSCRIPTS)
