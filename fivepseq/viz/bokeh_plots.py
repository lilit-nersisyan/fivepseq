import logging

import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import bokeh
from bokeh.colors import RGB
from bokeh.io import output_file, save, export_svgs
from bokeh.io import show
from bokeh.layouts import gridplot, row, widgetbox
from bokeh.models import HoverTool, Arrow, NormalHead, ColumnDataSource, LinearColorMapper, FactorRange, TableColumn, \
    DataTable, PanTool, BoxZoomTool, WheelZoomTool, SaveTool, ResetTool
from bokeh.plotting import figure
from bokeh.transform import transform
from fivepseq.logic.algorithms.count_stats.count_stats import CountStats

from fivepseq import config
from fivepseq.logic.structures import codons
from fivepseq.logic.structures.fivepseq_counts import CountManager

tools = [PanTool(), BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()]

def bokeh_composite(title, figure_list, filename, ncols=2):
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Making composite plot")

    output_file(title + ".html")

    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Title provided as: %s" % title)
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Number of figures: %d" % len(figure_list))
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Number of columns: %d" % ncols)

    p = gridplot(figure_list, ncols=ncols)
    print filename
    save(p, filename=filename)
    #p.output_backend = "svg"
    #export_svgs(p, title + ".svg")


def bokeh_table(title, table_df_dict):
    output_file(title + ".html")
    mainLayout = row(row(), name=title)
    for key in table_df_dict.keys():
        table_df = table_df_dict.get(key)
        table_dict = {}
        if table_df is not None:
            for i in range(table_df.shape[0]):
                table_dict.update({key + "_" + table_df.index[i]: list(table_df.iloc[i,:])})

            data_table = DataTable(source = ColumnDataSource(table_dict))
            p = figure(title = title + "_" + key, height = 500)

            mainLayout.children[0].children.append(widgetbox(data_table))
        else:
            mainLayout.children[0].children.append(None)

    return mainLayout

def bokeh_scatter_plot(title, region, count_series_dict, color_dict):
    if count_series_dict is None:
        return None
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Making count scatter plot: " + region)
    output_file(title + ".png")

    p = figure(title=title,
               x_axis_label="position from %s" % region, y_axis_label="5'seq read counts")

    for key in count_series_dict.keys():
        count_series = count_series_dict.get(key)
        if count_series is None:
            return None
        c = color_dict.get(key)
        line = p.line(count_series.D, count_series.C, legend=key, line_color=RGB(c[0], c[1], c[2]), line_width=2)

        hover = HoverTool(tooltips=[('position', '@x'), ('count', '@y')], renderers=[line])
        p.add_tools(hover)
    p.legend.click_policy = "hide"

    return p


def bokeh_fft_plot(title, align_region, signal_series_dict, color_dict, period_max=50):
    if signal_series_dict is None:
        return None
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Making FFT signal scatter plot: " + align_region)
    output_file(title + ".html")

    p = figure(title=title, x_range=(0, period_max),
               x_axis_label="Periodicity", y_axis_label="FFT signal")

    for key in signal_series_dict.keys():
        count_series = signal_series_dict.get(key)
        if count_series is None:
            return None
        c = color_dict.get(key)
        line = p.line(count_series.D, count_series.C, legend=key, line_color=RGB(c[0], c[1], c[2]),
                      line_width=2)
        hover = HoverTool(tooltips=[('period', '@x'), ('signal', '@y')], renderers=[line])
        p.add_tools(hover)
    p.legend.click_policy = "hide"
    return p


def bokeh_normalized_meta_scatter_plot(title, transcript_count_list_dict, color_dict,
                                       x_max, span_size, show_plot=False):
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Plotting normalized full-length meta counts. ")
    output_file(title + ".html", mode="cdn")

    p = figure(title=title,
               x_axis_label="normalized position within the transcript and surrounding span regions",
               y_axis_label="5'seq read counts")

    for key in transcript_count_list_dict.keys():
        count_vector_list = transcript_count_list_dict.get(key)
        # normalized_vector_list = [[]]*x_max
        normalized_count_sums = np.zeros(x_max)
        for index in range(0, len(count_vector_list)):
            count_vector = count_vector_list[index]
            stretch_vector = np.zeros(x_max * len(count_vector))
            for i in range(0, len(count_vector)):
                stretch_vector[(0 + i * x_max):((i + 1) * x_max)] = count_vector[i]
            normalized_vector = np.zeros(x_max)
            for i in range(0, x_max):
                normalized_vector[i] = np.mean(stretch_vector[(i * len(count_vector)):((i + 1) * len(count_vector))])
            # normalized_vector_list[index] = list(normalized_vector)
            normalized_count_sums = [sum(x) for x in zip(normalized_count_sums, normalized_vector)]
        normalized_meta_counts = [c / float(len(count_vector_list)) for c in normalized_count_sums]

        c = color_dict.get(key)
        line = p.line(list(range(0, x_max)), normalized_meta_counts, legend=key, line_color=RGB(c[0], c[1], c[2]))
        hover = HoverTool(tooltips=[('position', '@x'), ('count', '@y')], renderers=[line])
        p.add_tools(hover)
    p.legend.click_policy = "hide"

    if show_plot:
        show(p)
    return p


def bokeh_transcript_scatter_plot(title, transcript_count_list_dict, transcript_assembly, color_dict,
                                  align_to, span_size, index_filter, min_count=0, max_count=1000, save_plot=True):
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info(
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
                transcript = transcript_assembly.ID[index].split("transcript:")[1]
                gene = transcript_assembly.gene[index].split("gene:")[1]

                c = color_dict.get(key)
                line = p.line(count_series.index, count_series.values, legend=key, line_color=RGB(c[0], c[1], c[2]))

                p.add_tools(HoverTool(tooltips=[('position', '@x'), ('count', '@y'),
                                                ['transcript', transcript],
                                                ['gene', gene]], renderers=[line]))
    p.legend.click_policy = "hide"
    if save_plot:
        show(p)
    return p


def bokeh_triangle_plot(title, frame_df_dict, color_dict, transcript_index=None):
    if color_dict is None:
        return None
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Plotting triangle plots")
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
            [a, b, c] = frame_df.iloc[point, :][['F0', 'F1', 'F2']]
            if a == b == c:
                pass
            else:
                x[counter] = triangle_transform(a, b, c)[0]
                y[counter] = triangle_transform(a, b, c)[1]
                counter = counter + 1
        x = x[0:(counter - 1)]
        y = y[0:(counter - 1)]
        p.circle(x, y, color=color_dict.get(key), legend=key)

    p.legend.click_policy = "hide"
    return p


def triangle_transform(a, b, c):
    if a + b + c == 0:
        return -1, -1

    x = 1 / 2. * (float(a + 2. * b) / (a + b + c))
    y = (np.sqrt(3) / 2.) * (float(a) / (a + b + c))

    return x, y


def bokeh_heatmap_grid(title_prefix, amino_acid_df_dict, scale=False):
    if amino_acid_df_dict is None:
        return None

    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Plotting amino acid pauses: %s" % title_prefix)
    mainLayout = row(row(), name=title_prefix + ' amino acid pauses')

    for key in amino_acid_df_dict.keys():
        logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info(key)
        amino_acid_df = amino_acid_df_dict.get(key)
        if amino_acid_df is not None:
            if scale:
                for i in range(amino_acid_df.shape[0]):
                    amino_acid_df.iloc[i, :] /= sum(amino_acid_df.iloc[i, :])

            colormap = cm.get_cmap("viridis")
            bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
            mapper = LinearColorMapper(palette=bokehpalette, low=0, high=amino_acid_df.max().max())

            amino_acid_df.index.name = "aa"
            amino_acid_df.columns.name = "dist"
            df = amino_acid_df.stack().rename("value").reset_index()
            source = ColumnDataSource(df)

            p = figure(title=title_prefix + "_" + key,
                       x_range=FactorRange(factors=list(amino_acid_df.columns)),
                       y_range=FactorRange(factors=list(amino_acid_df.index)),
                       x_axis_label="distance from amino acid", y_axis_label="5'seq read counts")

            rect = p.rect(x='dist', y='aa', width=1, height=1,
                          source=source, fill_color=transform('value', mapper),
                          line_color=None)
            hover = HoverTool(tooltips=[('distance', '@dist'), ('count', '@value')], renderers=[rect])
            p.add_tools(hover)

            mainLayout.children[0].children.append(p)

    return mainLayout


def bokeh_frame_barplots(title_prefix, frame_df_dict, frame_stats_df_dict, color_dict):
    if frame_df_dict is None:
        return None
    logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).info("Plotting frame barplots")
    mainLayout = row(row(), name=title_prefix + ' frame histograms')

    counts_dict = {}
    count_max = 0
    for key in frame_df_dict.keys():

        frame_stats_df = frame_stats_df_dict.get(key)
        logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).debug("key: %s\n%s" % (key, str(frame_stats_df)))
        if frame_stats_df is not None:

            counts = [frame_stats_df.loc[CountStats.FRAME_COUNT, CountStats.F0],
                      frame_stats_df.loc[CountStats.FRAME_COUNT, CountStats.F1],
                      frame_stats_df.loc[CountStats.FRAME_COUNT, CountStats.F2]]

            if count_max < max(counts):
                count_max = max(counts)
            counts_dict.update({key: counts})
        else:
            counts_dict.update({key: None})

    for key in frame_df_dict.keys():
        color = color_dict.get(key)
        counts = counts_dict.get(key)
        if counts is None:
            #mainLayout.children[0].children.append(None)
            logging.getLogger(config.FIVEPSEQ_PLOT_LOGGER).warn("Frame counts stats not found for sample %s" % key)
        else:
            frames = ["F1", "F2", "F3"]
            p = figure(x_range=frames, y_range=(0, count_max),
                       plot_height=500, title=key + "_" + title_prefix)
            legend = ""
            frame_stats_df = frame_stats_df_dict.get(key)
            if frame_stats_df is None:
                legend = "Index: %.3f" % (np.log(float(counts[1] / np.mean([counts[0], counts[2]]))))
            else:
                legend = "p-vals:\tF0: %.2f\tF1: %.2f\tF2: %.2f" % (
                    np.round(frame_stats_df.loc[CountStats.PVAL_PAIR_MAX, CountStats.F0], 2),
                    np.round(frame_stats_df.loc[CountStats.PVAL_PAIR_MAX, CountStats.F1], 2),
                    np.round(frame_stats_df.loc[CountStats.PVAL_PAIR_MAX, CountStats.F2], 2),
                )
                p.vbar(x=frames, width=0.8, top=counts, bottom=0,
                       fill_color=color, line_color=None,
                       legend=legend)
                mainLayout.children[0].children.append(p)

    if len(mainLayout.children[0].children) == 0:
        return None
    return mainLayout


def bokeh_aa_scatter_grid(title_prefix, amino_acid_df_dict):
    print("Plotting amino acid pauses: %s" % title_prefix)
    # output_file(title_prefix + ".html")
    mainLayout = row(row(), name=title_prefix + ' amino acid pauses')

    for key in amino_acid_df_dict.keys():
        print(key)
        amino_acid_df = amino_acid_df_dict.get(key)
        p = figure(title=title_prefix + "_" + key,
                   x_axis_label="distance from amino acid", y_axis_label="5'seq read counts")

        for aa in amino_acid_df.index:
            c = codons.Codons.AMINO_ACID_COLORS.get(aa)
            line = p.line(amino_acid_df.columns, amino_acid_df.loc[aa,], legend=aa,
                          line_color=c)
            print line
            hover = HoverTool(tooltips=[('distance', '@x'), ('count', '@y')], renderers=[line])
            p.add_tools(hover)

        p.legend.click_policy = "hide"
        mainLayout.children[0].children.append(p)

    return mainLayout


def bokeh_aa_pause_scatterplot(title, amino_acid_df):
    p = figure(title=title,
               x_axis_label="distance from amino acid", y_axis_label="5'seq read counts")

    for aa in amino_acid_df.index:
        c = codons.Codons.AMINO_ACID_COLORS.get(aa)
        line = p.line(amino_acid_df.columns, amino_acid_df.loc(aa, ), legend=aa,
                      line_color=RGB(c[0], c[1], c[2]))
        hover = HoverTool(tooltips=[('distance', '@x'), ('count', '@y')], renderers=[line])
        p.add_tools(hover)

    p.legend.click_policy = "hide"
    return p