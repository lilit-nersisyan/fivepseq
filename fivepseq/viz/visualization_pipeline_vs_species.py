import os

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

group = "vs"

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
    #"Escherichia_coli",
    "Gardnerella_vaginalis",
    "Lactobacillus_crispatus",
    #"Lactobacillus_jensenii",
    #"Prevotella_fusca",
    #"Ureaplasma_parvum",

)
for org in organisms:
    dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_fivepseq_" + org
    dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + group + "/" + group + "_fivepseq_" + org + "_plots"

    if not os.path.exists(dir_5pseq_plots):
        os.mkdir(dir_5pseq_plots)

    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
