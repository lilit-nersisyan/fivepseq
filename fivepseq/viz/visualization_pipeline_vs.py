import os

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "vs"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_best_fivepseq_plots"

if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "VS_106_S25",
    "VS_153_S26",
    "VS_156_S27",
    "VS_162_S28",
    "VS_189_S29",
    "VS_35_S22",
    "VS_64_S23",
    "VS_88_S24",
)

organisms = (
    "Alteromonas_australica",
    "Bacillus_subtilis",
    "Escherichia_coli",
    "Gardnerella_vaginalis",
    "Lactobacillus_amylophilus",
    "Lactobacillus_crispatus",
    "Lactobacillus_jensenii",
    "Mageeibacillus_indolicus",
    "Prevotella_denticola",
    "Prevotella_fusca",
    "Prevotella_sp_oral_taxon_299",
    "Ruminococcus_champanellensis",
    "Salmonella_enterica",
    "Serratia_marcescens",
    "Ureaplasma_parvum",
    "Xanthomonas_euvesicatoria",

)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
