import os
from sets import Set

import pandas as pd
import numpy as np
import colorlover as cl
from bokeh.layouts import gridplot
from bokeh.models import HoverTool, Arrow, NormalHead, ColumnDataSource

from bokeh.plotting import figure
from bokeh.io import output_file
from bokeh.colors import RGB
from bokeh.io import show
from fivepseq import config
from fivepseq.logic.structures.counts import FivePSeqCounts
from fivepseq.logic.structures.fivepseq_counts import CountManager


def bokeh_composite(title, figure_list, ncols=2):
    output_file(title + ".html")

    p = gridplot(figure_list, ncols=ncols)
    show(p)


def bokeh_scatter_plot(title, region, count_dict, color_dict):
    output_file(title + ".html")

    p = figure(title=title,
               x_axis_label="position from %s" % region, y_axis_label="5'seq read counts")

    for key in count_dict.keys():
        count_series = count_dict.get(key)
        c = color_dict.get(key)
        line = p.line(count_series.D, count_series.C, legend=key, line_color=RGB(c[0], c[1], c[2]))
        print line
        hover = HoverTool(tooltips=[('position', '@x'), ('count', '@y')], renderers=[line])
        p.add_tools(hover)
    return p


def bokeh_transcript_scatter_plot(title, transcript_count_list_dict, transcript_assembly_dict, color_dict,
                                  align_to, span_size, index_filter, min_count=0, max_count=1000):
    config.logger.debug(
        "Plotting transcript specific counts. %d filtered transcript indices specified" % len(index_filter))
    output_file(title + "_minCount-" + str(min_count) + "_maxCount-" + str(max_count) + ".html", mode="cdn")

    p = figure(title=title,
               y_range=(0, 200),
               x_axis_label="position from %s" % align_to, y_axis_label="5'seq read counts")
    # try setting range limits - p = figure(x_range=[0, 10], y_range=(10, 20))

    for key in transcript_count_list_dict.keys():
        count_vector_list = transcript_count_list_dict.get(key)
        if index_filter is None:
            index_filter = range(0, len(count_vector_list))
        for index in index_filter:
            count_vector = count_vector_list[index]
            if min_count < sum(count_vector) < max_count:
                count_series = CountManager.count_vector_to_series(count_vector, align_to, tail=span_size)
                transcript = transcript_assembly_dict.get(key).ID[index].split("transcript:")[1]
                gene = transcript_assembly_dict.get(key).gene[index].split("gene:")[1]

                c = color_dict.get(key)
                line = p.line(count_series.index, count_series.values, legend=key, line_color=RGB(c[0], c[1], c[2]))

                p.add_tools(HoverTool(tooltips=[('position', '@x'), ('count', '@y'),
                                                ['transcript', transcript],
                                                ['gene', gene]], renderers=[line]))
    return p


def bokeh_triangle_plot(title, frame_df_dict, color_dict, transcript_index=None):
    print("Plotting triangle plots")
    if transcript_index is not None:
        print("%d filtered indices specified " % len(transcript_index))

    output_file(title + ".html")

    p = figure(title=title)

    # draw the axes
    f0 = (0, 0, 1)
    f1 = (0, 1, 0)
    f2 = (1, 0, 0)

    p.add_layout(Arrow(end=NormalHead(fill_color="gray"),
                       x_start=triangle_transform(*f0)[0], y_start=triangle_transform(*f0)[1],
                       x_end=triangle_transform(*f1)[0], y_end=triangle_transform(*f1)[1]))  # F0
    p.add_layout(Arrow(end=NormalHead(fill_color="gray"),
                       x_start=triangle_transform(*f1)[0], y_start=triangle_transform(*f1)[1],
                       x_end=triangle_transform(*f2)[0], y_end=triangle_transform(*f2)[1]))  # F1
    p.add_layout(Arrow(end=NormalHead(fill_color="gray"),
                       x_start=triangle_transform(*f2)[0], y_start=triangle_transform(*f2)[1],
                       x_end=triangle_transform(*f0)[0], y_end=triangle_transform(*f0)[1]))  # F2

    # draw points
    for key in frame_df_dict.keys():
        frame_df = frame_df_dict.get(key)
        if transcript_index is None:
            transcript_index = range(0, len(frame_df))
        frame_df = frame_df.iloc[transcript_index, :]
        x = [0] * len(frame_df)
        y = [0] * len(frame_df)
        counter = 0
        for point in range(0, len(frame_df)):
            [a, b, c] = frame_df.iloc[point, :][['F0','F1','F2']]
            if a == b == c:
                pass
            else:
                x[counter] = triangle_transform(a, b, c)[0]
                y[counter] = triangle_transform(a, b, c)[1]
                counter = counter + 1
        x = x[0:(counter - 1)]
        y = y[0:(counter - 1)]
        p.circle(x, y, color=color_dict.get(key), legend=key)
    return p


def triangle_transform(a, b, c):
    if a + b + c == 0:
        return -1, -1

    x = 1 / 2. * (float(a + 2. * b) / (a + b + c))
    y = (np.sqrt(3) / 2.) * (float(a) / (a + b + c))

    return x, y


########################################
#               DATA
########################################

dir_5pseq_microbiome = "/proj/sllstore2017018/lilit/5pseq_microbiome"
meta_count_dict_lpla = {
    "lpla_untreated": pd.read_csv(
        os.path.join(dir_5pseq_microbiome, "fivepseq_lpla_untreated", "meta_counts_TERM.txt"),
        sep="\t"),
    "lpla_fragmented": pd.read_csv(
        os.path.join(dir_5pseq_microbiome, "fivepseq_l_pla_fragmented", "meta_counts_TERM.txt"),
        sep="\t")
}
colors_dict_lpla = dict(zip(meta_count_dict_lpla.keys(),
                            np.array(cl.scales['3']['div']['RdGy'])[[0, 2]]))

dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"

meta_count_Hela_rep1_term_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_TERM.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_TERM.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_TERM.txt"),
        sep="\t", header=None, names=["D", "C"]),
}

meta_count_Hela_rep1_start_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_START.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_START.txt"),
        sep="\t", header=None, names=["D", "C"]),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_START.txt"),
        sep="\t", header=None, names=["D", "C"]),
}

colors_dict_Hela = dict(zip(meta_count_Hela_rep1_term_dict.keys(), cl.to_numeric(cl.scales['6']['qual']['Set1'])))

frame_count_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "frame_counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "frame_counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "frame_counts_TERM.txt"),
        sep="\t"),
}

transcript_count_term_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_TERM.txt"),
        sep="\t"),
}

transcript_count_start_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_START.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_START.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_START.txt"),
        sep="\t"),
}

transcript_count_full_Hela_rep1_dict = {
    "Hela-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "counts_FULL_LENGTH.txt")),
    "HelaCHX-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_FULL_LENGTH.txt")),
    "HelaFrag-rep1": CountManager.read_counts_as_list(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "counts_FULL_LENGTH.txt"))
}

transcript_assembly_Hela_rep1_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "transcript_assembly.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "transcript_assembly.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "transcript_assembly.txt"),
        sep="\t"),
}

#########################################
#          filter transcripts           #
#########################################

transcript_assembly = pd.read_csv(
    os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "transcript_assembly.txt"),
    sep="\t")

filtered_index = []
gene_set = Set()
unique_transcripts = []
for i in range(0, len(transcript_assembly)):
    gene = transcript_assembly.gene[i].split("gene:")[1]
    gene_set.add(gene)

gene_transcript_dict = {}
for i in range(0, len(transcript_assembly)):
    gene = transcript_assembly.gene[i].split("gene:")[1]
    transcript = transcript_assembly.ID[i].split("transcript:")[1]
    if gene_transcript_dict.has_key(gene):
        gene_transcript_dict.get(gene).append(transcript)
    else:
        gene_transcript_dict.update({gene: [transcript]})

for i in range(0, len(transcript_assembly)):
    gene = transcript_assembly.gene[i].split("gene:")[1]
    if len(gene_transcript_dict.get(gene)) == 1:
        unique_transcripts.append(i)

p_scatter_term = bokeh_scatter_plot("Hela-rep1_term", FivePSeqCounts.TERM, meta_count_Hela_rep1_term_dict,
                                    colors_dict_Hela)
p_scatter_start = bokeh_scatter_plot("Hela-rep1_start", FivePSeqCounts.START, meta_count_Hela_rep1_start_dict,
                                     colors_dict_Hela)
p_triangle = bokeh_triangle_plot("Hela-rep1_triangle", frame_count_Hela_rep1_dict, colors_dict_Hela, transcript_index=unique_transcripts)
p_transcripts = bokeh_transcript_scatter_plot("Hela-rep1_transcript_counts",
                                              transcript_count_full_Hela_rep1_dict,
                                              transcript_assembly_Hela_rep1_dict,
                                              colors_dict_Hela,
                                              FivePSeqCounts.TERM, 500, unique_transcripts, 300, 500)
bokeh_composite("Hela-rep1", [p_scatter_start, p_scatter_term, p_triangle, p_transcripts], 3)
