import os
from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "gs"
org = "Bacillus_velezensis"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_fivepseq_" + org
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_fivepseq_" + org + "_plots"

if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "2GS_X2_Antrum_TOT_S21",
    "2GS_X4_Fundus_P16_S11",
    "GS_X3_Antrum_P16_S8",
    "GS_X3_MID_P16_S5",
    "2GS_X2_Fundus_TOT_S19",
    "2GS_X4_Fundus_P8_S10",
    "GS_X3_Antrum_P8_S7",
    "GS_X3_MID_P8_S4",
    "2GS_X2_MID_TOT_S20",
    "2GS_X4_Fundus_S16_S12",
    "GS_X3_Antrum_S16_S9",
    "GS_X3_MID_S16_S6",
    "2GS_X4_Antrum_P16_S17",
    "2GS_X4_MID_P16_S14",
    "GS_X3_Fundus_P16_S2",
    "2GS_X4_Antrum_P8_S16",
    "2GS_X4_MID_P8_S13",
    "GS_X3_Fundus_P8_S1",
    "2GS_X4_Antrum_S16_S18",
    "2GS_X4_MID_S16_S15",
    "GS_X3_Fundus_S16_S3"
)

organisms = (
    "Bacillus_velezensis",
)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
