import os
from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "ts"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"

if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "TS-1_CAM-1_S34_R1_001",
    "TS-1_Ctr-1_S32_R1_001",
    "TS-1_H2O-1_S35_R1_001",
    "TS-1_MUP-1_S33_R1_001",
    "TS-2_CAM-1_S38_R1_001",
    "TS-2_Ctr-1_S36_R1_001",
    "TS-2_H2O-1_S39_R1_001",
    "TS-2_MUP-1_S37_R1_001",
    "TS-3_CAM-1_S42_R1_001",
    "TS-3_Ctr-1_S40_R1_001",
    "TS-3_H2O-1_S43_R1_001",
    "TS-3_MUP-1_S41_R1_001"
)


organisms = (
    "Acidisphaera_rubrifaciens_HS-AP3",
    "Bacillus_cereus",
    "Bacillus_subtilis",
    "Bacillus_thuringiensis",
    "Bacillus_velezensis",
    "Candidatus_Solibacter_usitatus",
    "Clostridium_botulinum",
    "Conexibacter_woesei",
    "Escherichia_coli",
    "Gemmata",
    "Klebsiella_pneumoniae",
    "Lactobacillus_crispatus",
    "Salmonella_enterica",
    "Singulisphaera_acidiphila",
    "Thermodesulfovibrio_yellowstonii"
)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
