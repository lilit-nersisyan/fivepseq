import os

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "topo"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "TOPO10-pEC622CAM1",
    "TOPO10-pEC622CAM2",
    "TOPO10-pEC622Ctr1",
    "TOPO10-pEC622Ctr2"
)

organisms = (
    "Escherichia_coli",
)

for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
