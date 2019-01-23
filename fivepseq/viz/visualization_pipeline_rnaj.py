import os
from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "rnaj"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)



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
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots, "rnaj")
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots, "WT")
