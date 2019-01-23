import os

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "caulo"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "Caulo_cres-Ctr1_S63_R1_001",
    "Caulo_cres-Ctr3_S64_R1_001",
    "Caulo_cres-S1_S65_R1_001",
    "Caulo_cres-S3_S66_R1_001"
)


organisms = (
    "Caulobacter_vibrioides",
)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
