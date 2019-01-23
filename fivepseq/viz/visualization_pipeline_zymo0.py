import os

from fivepseq.viz.visualization_pipeline import plot_samples_for_organism

parent_group = "zymo"
group = "zymo0"
dir_5pseq_group = "/proj/sllstore2017018/lilit/vagaga2/" + parent_group + "/" + group + "_best_fivepseq"
dir_5pseq_plots = "/proj/sllstore2017018/lilit/vagaga2/" + parent_group + "/" + group + "_best_fivepseq_plots"
if not os.path.exists(dir_5pseq_plots):
    os.mkdir(dir_5pseq_plots)

samples = (
    "ZymoBiomicsD6300-MCS1",
    "ZymoBiomicsD6300-MCS2",
)

organisms = (
    "Bacillus_subtilis",
    "Cryptococcus_neoformans",
    "Enterococcus_faecalis",
    "Escherichia_coli",
    "Lactobacillus_fermentum",
    "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa",
    "Saccharomyces_cerevisiae",
    "Salmonella_enterica",
    "Staphylococcus_aureus",
    "Bacillus_subtilis",
    "Cryptococcus_neoformans",
    "Enterococcus_faecalis",
    "Escherichia_coli",
    "Lactobacillus_fermentum",
    "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa",
    "Saccharomyces_cerevisiae",
    "Salmonella_enterica",
    "Staphylococcus_aureus"
)
for org in organisms:
    plot_samples_for_organism(group, samples, org, dir_5pseq_group, dir_5pseq_plots)
