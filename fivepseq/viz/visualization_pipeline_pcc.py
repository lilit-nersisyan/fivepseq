import os

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "pcc"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "PCC6083-Ctr1-2_S67_R1_001",
    "PCC6083-HS1-2_S68_R1_001"
)

organisms = (
    "Synechocystis_PCC_6803",
)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
