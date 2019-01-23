import os
from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "lreu"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "lreu_10000",
    "lreu_1000",
)

organisms = (
    "Lactobacillus_plantarum",
    "Lactobacillus_reuteri"
)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
