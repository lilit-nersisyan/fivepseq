from __future__ import print_function
from __future__ import print_function

import logging
import os
import random

import colorlover as cl
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from bokeh.colors import RGB
from bokeh.io import output_file, save, export_svgs, export_png
from bokeh.io import show
from bokeh.layouts import gridplot, row, widgetbox
from bokeh.models import HoverTool, Arrow, NormalHead, ColumnDataSource, LinearColorMapper, FactorRange, DataTable, \
    PanTool, BoxZoomTool, WheelZoomTool, SaveTool, ResetTool, Div, Legend, Panel, Tabs, LabelSet
from bokeh.plotting import figure
from bokeh.transform import transform
# from sklearn.decomposition import PCA

from fivepseq import config
from fivepseq.logic.algorithms.count_stats.count_stats import CountStats
from fivepseq.logic.structures import codons
from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts
from fivepseq.viz.header_html import get_div_logo, get_div_footer, get_div_title

tools = [PanTool(), BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()]

LEGEND_POSITION = 'above'
LEGEND_ALIGNMENT = 'bottom_left'
COMBINED = "combined"


def bokeh_composite(title, figure_list, filename, ncols=2):
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making composite plot")

    output_file(title + ".html", title=title)

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Title provided as: %s" % title)
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Number of figures: %d" % len(figure_list))
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Number of columns: %d" % ncols)

    # generate rows and then merge into a grid so that figures are not scattered to the right side of the page

    i = 0
    rows = []
    for r in range(int(len(figure_list)/ncols)):
        fig_row = row(row(), name=title)
        for c in range(ncols):
            if i == len(figure_list):
                break
            fig = figure_list[i]
            if fig is not None:
                fig_row.children[0].children.append(fig)
            else:
                print(f"None for figure {i} at r: {r} and c: {c}")
            i += 1
        rows.append(fig_row)
    p = gridplot(rows, ncols=1, toolbar_location="left")
    #p = gridplot(figure_list, ncols=ncols, toolbar_location="left") < fivepseq v1.3.1
    div_logo = Div(text=get_div_logo())
    div_title = Div(text=get_div_title(title))
    div_footer = Div(text=get_div_footer())
    save([div_logo, div_title, p, div_footer], filename=filename)


def bokeh_table(title, table_df_dict):
    output_file(title + ".html")
    mainLayout = row(row(), name=title)
    for key in table_df_dict.keys():
        table_df = table_df_dict.get(key)
        table_dict = {}
        if table_df is not None:
            for i in range(table_df.shape[0]):
                table_dict.update({key + "_" + table_df.index[i]: list(table_df.iloc[i, :])})

            data_table = DataTable(source=ColumnDataSource(table_dict))
            p = figure(title=title + "_" + key, height=500)

            mainLayout.children[0].children.append(widgetbox(data_table))
        else:
            mainLayout.children[0].children.append(None)

    return mainLayout


def bokeh_tabbed_line_chart(title, region, group_count_series_dict_dict, color_dict,
                            scale=False, lib_size_dict_dict=None,
                            combine_sum=False, combine_weighted=False,
                            combine_color=None,
                            png_dir=None, svg_dir=None):
    if group_count_series_dict_dict is None:
        return None

    #logging.getLogger(config.FIVEPSEQ_LOGGER).info(
    #    "Making count line charts with geneset-tabs: " + title + ": " + region)

    tab_list = []
    for group in group_count_series_dict_dict.keys():
        # gs_count_series_dict = gs_filter(count_series_dict, gs_transcriptID_dict, gs)
        # TODO filter by geneset
        count_series_dict_gs = group_count_series_dict_dict[group]

        if lib_size_dict_dict is not None:
            lib_size_dict = lib_size_dict_dict[group]
        else:
            lib_size_dict = None

        p_group = bokeh_line_chart(get_key_title(title, group), region, count_series_dict_gs, color_dict, scale=scale,
                                   combine_sum=combine_sum, combine_weighted=combine_weighted,
                                   combine_color=combine_color,
                                   lib_size_dict=lib_size_dict, png_dir=png_dir, svg_dir=svg_dir)
        tab_list.append(Panel(child=p_group, title=group))

    tabs = Tabs(tabs=tab_list)
    return tabs


def bokeh_line_chart(title, region, count_series_dict, color_dict, scale=False, lib_size_dict=None,
                     combine_sum=False,
                     combine_weighted=False, combine_color=None, png_dir=None, svg_dir=None):
    if count_series_dict is None:
        return None

    if combine_weighted and combine_sum:
        e_msg = "Exception plotting line chart %s: the options %s and %s cannot be true at the same time: " \
                % (title + ": " + region, "combine_sum", "combine_weighted")
        logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
        return None

    if combine_weighted:
        suffix = COMBINED + "-weighted"
        count_series_dict = {suffix: CountManager.combine_count_series(count_series_dict,
                                                                       lib_size_dict=lib_size_dict,
                                                                       scale=scale)}

    elif combine_sum:
        suffix = COMBINED + "-sum"
        count_series_dict = {suffix: CountManager.combine_count_series(count_series_dict)}

    if combine_weighted or combine_sum:
        if lib_size_dict is None:
            e_msg = "Exception plotting line chart %s: the option %s or %s cannot be supplied with lib_size_dict of None: " \
                    % (title + ": " + region, "combine_sum", "combine_weighted")
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
            return None

        lib_size_dict = {suffix: sum(lib_size_dict.values())}

        c = combine_color
        if c is None:
            c = get_random_color()
        color_dict = {suffix: c}

    #logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making count line chart %s with options: scaled(%s), "
    #                                               "combine_sum(%s), comine_weighted(%s): " % (title + ": " +
    #                                                                                           region,
    #                                                                                           str(scale),
    #                                                                                           str(combine_sum),
    #                                                                                           str(
    #                                                                                               combine_weighted)))

    output_file(title + ".png")

    my_x_label = "position from %s" % region
    if scale:
        my_y_label = "5'seq RPM"
    else:
        my_y_label = "5'seq raw counts"

    p = figure(title=title,
               x_axis_label=my_x_label, y_axis_label=my_y_label)

    p_png = None
    p_svg = None
    if png_dir is not None:
        p_png = figure(title=title,
                       x_axis_label=my_x_label, y_axis_label=my_y_label)
    if svg_dir is not None:
        p_svg = figure(title=title,
                       x_axis_label=my_x_label, y_axis_label=my_y_label)

    legend_items = []
    legend_items_png = []
    legend_items_svg = []

    for key in count_series_dict.keys():
        count_series = count_series_dict.get(key)

        if count_series is None:
            # return None
            continue

        key_title = get_key_title(title, key)
        p_key_png = None
        p_key_svg = None
        if png_dir is not None:
            p_key_png = figure(title=key_title, x_axis_label=my_x_label, y_axis_label=my_y_label)
        if svg_dir is not None:
            p_key_svg = figure(title=key_title, x_axis_label=my_x_label, y_axis_label=my_y_label)
        c = color_dict.get(key)
        if c is None:
            logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Color not set for sample %s" % key)
            c = get_random_color()
        try:
            if scale and not combine_weighted:
                if lib_size_dict is None:
                    logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Lib size not specified: local counts will "
                                                                      "be used instead")
                    lib_size = sum(count_series.C) + 1
                else:
                    lib_size = lib_size_dict.get(key)
                y = (10 ** 6) * count_series.C / lib_size
            else:
                y = count_series.C
            source = ColumnDataSource(data=dict(
                x=list(count_series.D),
                y=list(y),
                raw=list(count_series.C),
                name=[key] * len(y)
            ))
            line = p.line('x', 'y', line_color=c, line_width=2, source=source)
            legend_items.append((key, [line]))
        except Exception as e:
            print("Error at key %s, reason; %s" % (key, str(e)))
            return None
        if scale:
            hover = HoverTool(tooltips=[('name', '@name'), ('position', '@x'), ('RPM', '@y'), ('raw count', '@raw')],
                              renderers=[line])
        else:
            hover = HoverTool(tooltips=[('name', '@name'), ('position', '@x'), ('raw count', '@raw')], renderers=[line])
        p.add_tools(hover)

        # figures for exporting
        if p_key_png is not None:
            p_key_png.line('x', 'y', line_color=c, line_width=2, source=source)
            export_images(p_key_png, key_title, png_dir=png_dir)
        if p_key_svg is not None:
            p_key_svg.line('x', 'y', line_color=c, line_width=2, source=source)
            export_images(p_key_svg, key_title, svg_dir=svg_dir)

        if p_png is not None:
            line_png = p_png.line('x', 'y', line_color=c, line_width=2, source=source)
            legend_items_png.append((key, [line_png]))
        if p_svg is not None:
            line_svg = p_svg.line('x', 'y', line_color=c, line_width=2, source=source)
            legend_items_svg.append((key, [line_svg]))

    legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
    p.add_layout(legend, LEGEND_POSITION)
    p.legend.click_policy = "hide"

    if p_png is not None:
        legend_png = Legend(items=legend_items_png, location=LEGEND_ALIGNMENT)
        p_png.add_layout(legend_png, LEGEND_POSITION)
        export_images(p_png, title + "_overlay", png_dir=png_dir)

    if p_svg is not None:
        legend_svg = Legend(items=legend_items_svg, location=LEGEND_ALIGNMENT)
        p_svg.add_layout(legend_svg, LEGEND_POSITION)
        export_images(p_svg, title + "_overlay", svg_dir=svg_dir)

    return p


def bokeh_tabbed_fft_plot(title, align_region, group_signal_series_dict_dict, color_dict, period_max=50,
                          lib_size_dict_dict=None, combine_sum=False, combine_weighted=False, combine_color=None,
                          png_dir=None,
                          svg_dir=None):
    if group_signal_series_dict_dict is None:
        return None

    logging.getLogger(config.FIVEPSEQ_LOGGER).info(
        "Making fft plots with geneset-tabs: " + title + ": " + align_region)

    tab_list = []
    for group in group_signal_series_dict_dict.keys():
        signal_series_dict_gs = group_signal_series_dict_dict[group]

        if lib_size_dict_dict is not None:
            lib_size_dict = lib_size_dict_dict[group]
        else:
            lib_size_dict = None

        p_group = bokeh_fft_plot(get_key_title(title, group), align_region, signal_series_dict_gs, color_dict,
                                 period_max=period_max,
                                 lib_size_dict=lib_size_dict,
                                 combine_sum=combine_sum, combine_weighted=combine_weighted,
                                 combine_color=combine_color,
                                 png_dir=png_dir, svg_dir=svg_dir)
        tab_list.append(Panel(child=p_group, title=group))

    tabs = Tabs(tabs=tab_list)
    return tabs


def bokeh_fft_plot(title, align_region, signal_series_dict, color_dict, period_max=50,
                   lib_size_dict=None, combine_sum=False, combine_weighted=False, combine_color=None,
                   png_dir=None, svg_dir=None):
    if signal_series_dict is None:
        return None

    if combine_weighted and combine_sum:
        e_msg = "Exception plotting line chart %s: the options %s and %s cannot be true at the same time: " \
                % (title + ": " + align_region, "combine_sum", "combine_weighted")
        logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
        return None

    if combine_weighted:
        suffix = COMBINED + "-weighted"
        signal_series_dict = {suffix: CountManager.combine_count_series(signal_series_dict, lib_size_dict)}

    elif combine_sum:
        suffix = COMBINED + "-sum"
        signal_series_dict = {suffix: CountManager.combine_count_series(signal_series_dict)}

    if combine_weighted or combine_sum:
        if lib_size_dict is None:
            e_msg = "Exception plotting line chart %s: the option %s or %s cannot be supplied with lib_size_dict of None: " \
                    % (title + ": " + align_region, "combine_sum", "combine_weighted")
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
            return None

        c = combine_color
        if c is None:
            c = get_random_color()
        color_dict = {suffix: c}

    #logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making count line chart %s with options: "
    #                                               "combine_sum(%s), comine_weighted(%s): "
    #                                               % (title + ": " +
    #                                                  align_region,
    #                                                  str(combine_sum),
    #                                                  str(combine_weighted)))

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making FFT signal line chart: " + align_region)
    output_file(title + ".html")

    my_x_label = "Periodicity"
    my_y_label = "FFT signal"

    p = figure(title=title, x_range=(0, period_max), x_axis_label=my_x_label, y_axis_label=my_y_label)

    p_png = None
    p_svg = None
    if png_dir is not None:
        p_png = figure(title=title, x_range=(0, period_max), x_axis_label=my_x_label, y_axis_label=my_y_label)
    if svg_dir is not None:
        p_svg = figure(title=title, x_range=(0, period_max), x_axis_label=my_x_label, y_axis_label=my_y_label)

    legend_items = []
    legend_items_png = []
    legend_items_svg = []
    for key in signal_series_dict.keys():
        count_series = signal_series_dict.get(key)
        if count_series is None:
            return None
        key_title = get_key_title(title, key)
        c = color_dict.get(key)
        if c is None:
            logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Color not set for sample %s" % key)
            c = get_random_color()

        source = ColumnDataSource(data=dict(
            x=list(count_series.D),
            y=list(count_series.C),
            name=[key] * len(count_series.C),
        ))

        line = p.line('x', 'y', line_color=c, line_width=2, source=source)
        hover = HoverTool(tooltips=[('name', '@name'), ('period', '@x'), ('signal', '@y')], renderers=[line])
        p.add_tools(hover)
        legend_items.append((key, [line]))

        # figures for exporting
        if png_dir is not None:
            p_key_png = figure(title=key_title, x_axis_label=my_x_label, y_axis_label=my_y_label)
            p_key_png.line(count_series.D, count_series.C, line_color=c, line_width=2)
            export_images(p_key_png, key_title, png_dir=png_dir)
        if svg_dir is not None:
            p_key_svg = figure(title=key_title, x_axis_label=my_x_label, y_axis_label=my_y_label)
            p_key_svg.line(count_series.D, count_series.C, line_color=c, line_width=2)
            export_images(p_key_svg, key_title, svg_dir=svg_dir)

        if p_png is not None:
            line_png = p_png.line(count_series.D, count_series.C, line_color=c, line_width=2)
            legend_items_png.append((key, [line_png]))
        if p_svg is not None:
            line_svg = p_svg.line(count_series.D, count_series.C, line_color=c, line_width=2)
            legend_items_svg.append((key, [line_svg]))

    legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
    p.add_layout(legend, LEGEND_POSITION)
    p.legend.click_policy = "hide"

    if p_png is not None:
        legend_png = Legend(items=legend_items_png, location=LEGEND_ALIGNMENT)
        p_png.add_layout(legend_png, LEGEND_POSITION)
        export_images(p_png, title + "_overlay", png_dir=png_dir)

    if p_svg is not None:
        legend_svg = Legend(items=legend_items_svg, location=LEGEND_ALIGNMENT)
        p_svg.add_layout(legend_svg, LEGEND_POSITION)
        export_images(p_svg, title + "_overlay", svg_dir=svg_dir)

    return p


def bokeh_normalized_meta_line_chart(title, transcript_count_list_dict, color_dict,
                                     x_max, show_plot=False, png_dir=None, svg_dir=None):
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Plotting normalized full-length meta counts. ")
    output_file(title + ".html", mode="cdn")

    p = figure(title=title,
               x_axis_label="normalized position within the transcript and surrounding span regions",
               y_axis_label="5'seq read counts")

    legend_items = []
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

        if c is None:
            logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Color not set for sample %s" % key)
            c = get_random_color()

        line = p.line(list(range(0, x_max)), normalized_meta_counts, line_color=RGB(c[0], c[1], c[2]))
        legend_items.append((key, [line]))
        hover = HoverTool(tooltips=[('position', '@x'), ('count', '@y')], renderers=[line])
        p.add_tools(hover)

    legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
    p.add_layout(legend, LEGEND_POSITION)
    p.legend.click_policy = "hide"

    if show_plot:
        show(p)
    return p


def bokeh_transcript_line_chart(title, transcript_count_list_dict, transcript_assembly, color_dict,
                                align_to, span_size, index_filter, min_count=0, max_count=1000, save_plot=True):
    logging.getLogger(config.FIVEPSEQ_LOGGER).info(
        "Plotting transcript specific counts. %d filtered transcript indices specified" % len(index_filter))
    output_file(title + "_minCount-" + str(min_count) + "_maxCount-" + str(max_count) + ".html", mode="cdn")

    p = figure(title=title,
               y_range=(0, 200),
               x_axis_label="position from %s" % align_to, y_axis_label="5'seq read counts")
    # try setting range limits - p = figure(x_range=[0, 10], y_range=(10, 20))

    legend_items = []
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
                if c is None:
                    logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Color not set for sample %s" % key)
                    c = cl.scales['9']['qual']['Set3'][1]
                # line = p.line(count_series.index, count_series.values, line_color=RGB(c[0], c[1], c[2]))
                line = p.line(count_series.index, count_series.values, line_color=c)
                legend_items.append((key, [line]))
                p.add_tools(HoverTool(tooltips=[('position', '@x'), ('count', '@y'),
                                                ['transcript', transcript],
                                                ['gene', gene]], renderers=[line]))

    legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
    p.add_layout(legend, LEGEND_POSITION)
    p.legend.click_policy = "hide"

    if save_plot:
        show(p)
    return p


def bokeh_tabbed_triangle_plot(title, group_frame_df_dict_dict, color_dict,
                               transcript_names_dict,
                               lib_size_dict_dict=None,
                               combine_sum=False, combine_weighted=False, combine_color=None,
                               png_dir=None, svg_dir=None, count_threshold=10):
    if group_frame_df_dict_dict is None:
        return None

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making tabbed triangle plots: " + title)

    tab_list = []
    for group in group_frame_df_dict_dict.keys():
        # gs_count_series_dict = gs_filter(count_series_dict, gs_transcriptID_dict, gs)
        # TODO filter by geneset
        frame_df_dict = group_frame_df_dict_dict[group]

        if lib_size_dict_dict is not None:
            lib_size_dict = lib_size_dict_dict[group]
        else:
            lib_size_dict = None

        if group in transcript_names_dict:
            group_transcript_names_dict = transcript_names_dict[group]
        else:
            group_transcript_names_dict = transcript_names_dict

        p_group = bokeh_triangle_plot(get_key_title(title, group), frame_df_dict, color_dict,
                                      transcript_names_dict=group_transcript_names_dict,
                                      lib_size_dict=lib_size_dict,
                                      combine_sum=combine_sum, combine_weighted=combine_weighted,
                                      combine_color=combine_color,
                                      png_dir=png_dir, svg_dir=svg_dir,
                                      count_threshold=count_threshold)
        tab_list.append(Panel(child=p_group, title=group))
    tabs = Tabs(tabs=tab_list)
    return tabs


def bokeh_triangle_plot(title, frame_df_dict, color_dict, transcript_names_dict, lib_size_dict=None,
                        combine_sum=False, combine_weighted=False, combine_color=None,
                        transcript_index_filter=None, png_dir=None, svg_dir=None, count_threshold=10):
    if color_dict is None:
        return None
    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Plotting triangle plots")
    if transcript_index_filter is not None:
        print("%d filtered indices specified " % len(transcript_index_filter))

    output_file(title + ".html")

    if combine_weighted and combine_sum:
        e_msg = "Exception plotting line chart %s: the options %s and %s cannot be true at the same time: " \
                % (title + ": ", "combine_sum", "combine_weighted")
        logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
        return None

    if combine_weighted:
        suffix = COMBINED + "-w"
        frame_df_dict = {suffix: CountManager.combine_frame_counts(frame_df_dict, lib_size_dict)}
        transcript_names_dict = {suffix: transcript_names_dict[next(iter(transcript_names_dict))]}
        # NOTE this may be a problem when the lists do not match (which they should not) - but I'm lazy to check for this

    elif combine_sum:
        suffix = COMBINED + "-s"
        frame_df_dict = {suffix: CountManager.combine_frame_counts(frame_df_dict, lib_size_dict)}
        transcript_names_dict = {suffix: transcript_names_dict[next(iter(transcript_names_dict))]}

    if combine_weighted or combine_sum:
        if lib_size_dict is None:
            e_msg = "Exception plotting line chart %s: the option %s or %s cannot be supplied with lib_size_dict of None: " \
                    % (title, "combine_sum", "combine_weighted")
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
            return None

        c = combine_color
        if c is None:
            c = get_random_color()
        color_dict = {suffix: c}

    p = get_empty_triangle_canvas(title=title)

    # figures for export
    p_png = None
    p_svg = None
    if png_dir is not None:
        p_png = get_empty_triangle_canvas(title=title)
    if svg_dir is not None:
        p_svg = get_empty_triangle_canvas(title=title)

    # draw points
    legend_items = []
    legend_items_png = []
    legend_items_svg = []
    for key in frame_df_dict.keys():
        key_title = get_key_title(title, key)
        p_key_png = None
        p_key_svg = None

        if png_dir is not None:
            p_key_png = get_empty_triangle_canvas(key_title)
        if svg_dir is not None:
            p_key_svg = get_empty_triangle_canvas(key_title)

        transcript_names = transcript_names_dict.get(key)

        frame_df = frame_df_dict.get(key)
        if transcript_index_filter is None:
            transcript_index = range(0, len(frame_df))
        else:
            transcript_index = transcript_index_filter

        frame_df = frame_df.iloc[transcript_index, :]
        x = [0] * len(frame_df)
        y = [0] * len(frame_df)
        f0 = [0] * len(frame_df)
        f1 = [0] * len(frame_df)
        f2 = [0] * len(frame_df)
        sample = [key] * len(frame_df)
        name = [key] * len(frame_df)
        counter = 0
        for point in range(0, len(frame_df)):
            [a, b, c] = frame_df.iloc[point, :][['F0', 'F1', 'F2']]
            if a == b == c == 0:
                # TODO why am I passing here (added the == 0 at the end)?
                pass
            if sum([a, b, c]) < count_threshold:  # TODO makes this an option in the future
                pass
            else:
                x[counter] = triangle_transform(a, b, c)[0]
                y[counter] = triangle_transform(a, b, c)[1]
                f0[counter] = frame_df.loc[point, 'F0']
                f1[counter] = frame_df.loc[point, 'F1']
                f2[counter] = frame_df.loc[point, 'F2']
                name[counter] = transcript_names[point]
                counter = counter + 1
        x = x[0:counter]
        y = y[0:counter]
        f0 = f0[0:counter]
        f1 = f1[0:counter]
        f2 = f2[0:counter]
        sample = sample[0:counter]
        name = name[0:counter]

        source = ColumnDataSource(data=dict(
            x=x,
            y=y,
            F0=f0,
            F1=f1,
            F2=f2,
            name=name,
            sample=sample,
        ))

        # circles = p.circle(x, y, color=color_dict.get(key), fill_alpha=0.5, size=10, line_color=None)
        circles = p.circle('x', 'y', color=color_dict.get(key), fill_alpha=0.5, size=5, line_color=None, source=source)
        legend_items.append((key, [circles]))
        hover = HoverTool(
            tooltips=[('name', '@name'), ('sample', '@sample'), ('F0', '@F0'), ('F1', '@F1'), ('F2', '@F2')],
            renderers=[circles])
        p.add_tools(hover)

        if p_key_png is not None:
            circles_png = p_key_png.circle(x, y, color=color_dict.get(key), fill_alpha=0.5, size=5, line_color=None)
            legend_items_png.append((key, [circles_png]))
            export_images(p_key_png, key_title, png_dir=png_dir)
        if p_key_svg is not None:
            circles_svg = p_key_svg.circle(x, y, color=color_dict.get(key), fill_alpha=0.5, size=5, line_color=None)
            legend_items_svg.append((key, [circles_svg]))
            export_images(p_key_svg, key_title, svg_dir=svg_dir)

        if p_png is not None:
            p_png.circle(x, y, color=color_dict.get(key))
        if p_svg is not None:
            p_svg.circle(x, y, color=color_dict.get(key))

    legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
    p.add_layout(legend, LEGEND_POSITION)
    p.legend.click_policy = "hide"

    if p_png is not None:
        legend_png = Legend(items=legend_items_png, location=LEGEND_ALIGNMENT)
        p_png.add_layout(legend_png, LEGEND_POSITION)
        export_images(p_png, title + "_overlay", png_dir=png_dir)

    if p_svg is not None:
        legend_svg = Legend(items=legend_items_svg, location=LEGEND_ALIGNMENT)
        p_svg.add_layout(legend_svg, LEGEND_POSITION)
        export_images(p_svg, title + "_overlay", svg_dir=svg_dir)

    return p


def get_empty_triangle_canvas(title):
    p = figure(title=title, x_range=(-0.1, 1.1), y_range=(-0.1, 1.1))
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

    source = ColumnDataSource(
        data=dict(x=[triangle_transform(*f1)[0], triangle_transform(*f2)[0], triangle_transform(*f0)[0]],
                  y=[triangle_transform(*f1)[1], triangle_transform(*f2)[1], triangle_transform(*f0)[1]],
                  text=['F1', 'F0', 'F2'],
                  x_offset=[0, 0, -20]))
    p.add_layout(LabelSet(x='x', y='y', text='text',
                          x_offset='x_offset',
                          source=source, level='glyph', render_mode='canvas'))
    return p


def triangle_transform(a, b, c):
    if a + b + c == 0:
        return -1, -1

    x = 1 / 2. * (float(a + 2. * b) / (a + b + c))
    y = (np.sqrt(3) / 2.) * (float(a) / (a + b + c))

    return x, y


def bokeh_tabbed_heatmap_grid(title_prefix, group_amino_acid_df_dict_dict,
                              lib_size_dict_dict=None,
                              combine_sum=False, combine_weighted=False,
                              scale=False, png_dir=None, svg_dir=None):
    if group_amino_acid_df_dict_dict is None:
        return None

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making tabbed heatmap: " + title_prefix)

    tab_list = []
    for group in group_amino_acid_df_dict_dict.keys():
        amino_acid_df_dict = group_amino_acid_df_dict_dict[group]

        if lib_size_dict_dict is not None:
            lib_size_dict = lib_size_dict_dict[group]
        else:
            lib_size_dict = None

        p_group = bokeh_heatmap_grid(get_key_title(title_prefix, group), amino_acid_df_dict, scale=scale,
                                     lib_size_dict=lib_size_dict, combine_sum=combine_sum,
                                     combine_weighted=combine_weighted,
                                     png_dir=png_dir, svg_dir=svg_dir)
        tab_list.append(Panel(child=p_group, title=group))
    tabs = Tabs(tabs=tab_list)
    return tabs


def bokeh_heatmap_grid(title_prefix, amino_acid_df_dict, scale=False, lib_size_dict=None,
                       combine_sum=False, combine_weighted=False,
                       png_dir=None, svg_dir=None):
    if amino_acid_df_dict is None:
        return None

    if combine_weighted and combine_sum:
        e_msg = "Exception plotting heatmap %s: the options %s and %s cannot be true at the same time: " \
                % (title_prefix, "combine_sum", "combine_weighted")
        logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
        return None

    if combine_weighted:
        suffix = COMBINED + "-weighted"
        amino_acid_df_dict = {suffix: CountManager.combine_amino_acid_dfs(amino_acid_df_dict, lib_size_dict)}

    elif combine_sum:
        suffix = COMBINED + "-sum"
        amino_acid_df_dict = {suffix: CountManager.combine_amino_acid_dfs(amino_acid_df_dict)}

    if combine_weighted or combine_sum:
        if lib_size_dict is None:
            e_msg = "Exception plotting line chart %s: the option %s or %s cannot be supplied with lib_size_dict of None: " \
                    % (title_prefix, "combine_sum", "combine_weighted")
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
            return None

    #logging.getLogger(config.FIVEPSEQ_LOGGER).info(
    #    "Making heatmaps for amino acid relative counts for %s with options: scaled(%s), "
    #    "combine_sum(%s), comine_weighted(%s): " % (title_prefix,
    #                                                str(scale),
    #                                                str(combine_sum),
    #                                                str(combine_weighted)))

    mainLayout = row(row(), name=title_prefix + ' amino acid pauses')

    if not scale:
        y_axis_label = "5'seq raw read counts"
    else:
        y_axis_label = "5'seq scaled counts"

    for key in amino_acid_df_dict.keys():
        logging.getLogger(config.FIVEPSEQ_LOGGER).info(key)
        amino_acid_df = amino_acid_df_dict.get(key).copy(deep=True)
        if amino_acid_df is not None and len(amino_acid_df) > 0:
            if scale:
                for i in range(amino_acid_df.shape[0]):
                    amino_acid_df.iloc[i, :] /= sum(amino_acid_df.iloc[i, :]) + 1
            #                    amino_acid_df.iloc[i, :] -= amino_acid_df.iloc[i, :].min()
            #                    if amino_acid_df.iloc[i,:].max() > 0:
            #                        amino_acid_df.iloc[i, :] /= amino_acid_df.iloc[i,:].max()

            colormap = cm.get_cmap("viridis")
            bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
            mapper = LinearColorMapper(palette=bokehpalette, low=0, high=amino_acid_df.max().max())

            amino_acid_df.index.name = "aa"
            amino_acid_df.columns.name = "dist"
            df = amino_acid_df.stack().rename("value").reset_index()
            source = ColumnDataSource(df)

            key_title = get_key_title(title_prefix, key)
            p = figure(title=key_title,
                       x_range=FactorRange(factors=list(amino_acid_df.columns)),
                       y_range=FactorRange(factors=list(amino_acid_df.index)),
                       x_axis_label="distance from amino acid", y_axis_label=y_axis_label)

            # figures for export
            p_png = None
            p_svg = None
            if png_dir is not None:
                p_png = figure(title=key_title,
                               x_range=FactorRange(factors=list(amino_acid_df.columns)),
                               y_range=FactorRange(factors=list(amino_acid_df.index)),
                               x_axis_label="distance from amino acid", y_axis_label=y_axis_label)
            if svg_dir is not None:
                p_svg = figure(title=key_title,
                               x_range=FactorRange(factors=list(amino_acid_df.columns)),
                               y_range=FactorRange(factors=list(amino_acid_df.index)),
                               x_axis_label="distance from amino acid", y_axis_label=y_axis_label)

            rect = p.rect(x='dist', y='aa', width=1, height=1,
                          source=source, fill_color=transform('value', mapper),
                          line_color=None)
            hover = HoverTool(tooltips=[('distance', '@dist'), ('count', '@value')], renderers=[rect])
            p.add_tools(hover)

            if p_png is not None:
                p_png.rect(x='dist', y='aa', width=1, height=1,
                           source=source, fill_color=transform('value', mapper),
                           line_color=None)
            if p_svg is not None:
                p_svg.rect(x='dist', y='aa', width=1, height=1,
                           source=source, fill_color=transform('value', mapper),
                           line_color=None)

            mainLayout.children[0].children.append(p)
            export_images(p_png, key_title, png_dir=png_dir)
            export_images(p_svg, key_title, svg_dir=svg_dir)

    return mainLayout


def bokeh_tabbed_frame_barplots(title, group_frame_df_dict_dict, group_frame_stats_df_dict_dict, color_dict,
                                scale=False, lib_size_dict_dict=None, combine_sum=False,
                                combine_weighted=False, combine_color=None,
                                png_dir=None, svg_dir=None):
    if group_frame_df_dict_dict is None:
        return None

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making tabbed triangle plots: " + title)

    tab_list = []
    for group in group_frame_df_dict_dict.keys():
        # gs_count_series_dict = gs_filter(count_series_dict, gs_transcriptID_dict, gs)
        # TODO filter by geneset
        frame_df_dict = group_frame_df_dict_dict[group]
        if group_frame_stats_df_dict_dict is not None:
            frame_stats_df_dict = group_frame_stats_df_dict_dict[group]
        else:
            frame_stats_df_dict = {}
            for sample in frame_df_dict.keys():
                frame_stats_df_dict.update({sample: None})

        if lib_size_dict_dict is not None:
            lib_size_dict = lib_size_dict_dict[group]
        else:
            lib_size_dict = None

        p_group = bokeh_frame_barplots(get_key_title(title, group), frame_df_dict, frame_stats_df_dict, color_dict,
                                       lib_size_dict=lib_size_dict, scale=scale,
                                       combine_sum=combine_sum,
                                       combine_weighted=combine_weighted,
                                       png_dir=png_dir,
                                       svg_dir=svg_dir)
        tab_list.append(Panel(child=p_group, title=group))
    tabs = Tabs(tabs=tab_list)
    return tabs


def bokeh_frame_barplots(title_prefix, frame_df_dict, frame_stats_df_dict, color_dict,
                         lib_size_dict=None, scale=False,
                         combine_sum=False,
                         combine_weighted=False, combine_color=None,
                         png_dir=None, svg_dir=None):
    if frame_df_dict is None:
        return None

    if combine_weighted and combine_sum:
        e_msg = "Exception plotting line chart %s: the options %s and %s cannot be true at the same time: " \
                % (title_prefix, "combine_sum", "combine_weighted")
        logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
        return None

    if combine_weighted:
        suffix = COMBINED + "-weighted"
        frame_df_dict = {suffix: CountManager.combine_frame_counts(frame_df_dict, lib_size_dict)}

    elif combine_sum:
        suffix = COMBINED + "-sum"
        frame_df_dict = {suffix: CountManager.combine_frame_counts(frame_df_dict)}

    if combine_weighted or combine_sum:
        if lib_size_dict is None:
            e_msg = "Exception plotting line chart %s: the option %s or %s cannot be supplied with lib_size_dict of None: " \
                    % (title_prefix, "combine_sum", "combine_weighted")
            logging.getLogger(config.FIVEPSEQ_LOGGER).error(e_msg)
            return None

        c = combine_color
        if c is None:
            c = get_random_color()
        color_dict = {suffix: c}

    logging.getLogger(config.FIVEPSEQ_LOGGER).info("Making frame barplot %s with options: "
                                                   "combine_sum(%s), combine_weighted(%s): " % (title_prefix,
                                                                                                str(combine_sum),
                                                                                                str(combine_weighted)))

    mainLayout = row(row(), name=title_prefix + ' frame histograms')

    counts_dict = {}
    count_max = 0
    for key in frame_df_dict.keys():

        if frame_stats_df_dict is not None:
            frame_stats_df = frame_stats_df_dict.get(key)
        else:
            frame_stats_df = None
        #logging.getLogger(config.FIVEPSEQ_LOGGER).debug("key: %s\n%s" % (key, str(frame_stats_df)))
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
        legend_items = []
        color = color_dict.get(key)
        counts = counts_dict.get(key)[:]  # copying
        if counts is None:
            # mainLayout.children[0].children.append(None)
            logging.getLogger(config.FIVEPSEQ_LOGGER).warn("Frame counts stats not found for sample %s" % key)
        else:
            frames = ["F0", "F1", "F2"]
            if scale:
                counts /= sum(counts)
                counts *= 100
                max_y = 100
                y_label = "% of genome-wide total counts on the frame"
            else:
                max_y = count_max
                y_label = "Genome-wide total counts on the frame"
            key_title = get_key_title(title_prefix, key)
            p = figure(x_range=frames, y_range=(0, max_y),
                       plot_height=500, title=key_title,
                       y_axis_label=y_label)

            # figures for export
            p_png = None
            p_svg = None
            if png_dir is not None:
                p_png = figure(x_range=frames, y_range=(0, max_y),
                               plot_height=1000, plot_width=1000, title=title_prefix)
            if svg_dir is not None:
                p_svg = figure(x_range=frames, y_range=(0, max_y),
                               plot_height=500, title=title_prefix)

            bars = p.vbar(x=frames, width=0.8, top=counts, bottom=0,
                          fill_color=color, line_color=None)

            frame_stats_df = frame_stats_df_dict.get(key)
            if frame_stats_df is None:
                legend_text = "Index: %.3f" % (np.log(float(counts[1] / np.mean([counts[0], counts[2]]))))
                legend_items.append((legend_text, [bars]))
            else:

                legend_text = "FPI:    \tF0=%.2f | \tF1=%.2f | \tF2=%.2f" % (
                    np.round(frame_stats_df.loc[CountStats.FPI, CountStats.F0], 2),
                    np.round(frame_stats_df.loc[CountStats.FPI, CountStats.F1], 2),
                    np.round(frame_stats_df.loc[CountStats.FPI, CountStats.F2], 2),
                )
                legend_items.append((legend_text, [bars]))
                legend_text = "p-adj: \tF0=%.2f | \tF1=%.2f | \tF2=%.2f" % (
                    np.round(frame_stats_df.loc[CountStats.PVAL_PAIR_MAX, CountStats.F0], 2),
                    np.round(frame_stats_df.loc[CountStats.PVAL_PAIR_MAX, CountStats.F1], 2),
                    np.round(frame_stats_df.loc[CountStats.PVAL_PAIR_MAX, CountStats.F2], 2),
                )
                legend_items.append((legend_text, [bars]))

            legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
            p.add_layout(legend, LEGEND_POSITION)
            p.legend.click_policy = "hide"

            if p_png is not None:
                p_png.vbar(x=frames, width=0.8, top=counts, bottom=0,
                           fill_color=color, line_color=None)

            if p_svg is not None:
                p_png.vbar(x=frames, width=0.8, top=counts, bottom=0,
                           fill_color=color, line_color=None)

            mainLayout.children[0].children.append(p)
            export_images(p_png, key_title, png_dir=png_dir)
            export_images(p_svg, key_title, svg_dir=svg_dir)

    if len(mainLayout.children[0].children) == 0:
        return None
    return mainLayout


def bokeh_aa_line_grid(title_prefix, amino_acid_df_dict, png_dir=None, svg_dir=None):
    print("Plotting amino acid pauses: %s" % title_prefix)
    # output_file(title_prefix + ".html")
    mainLayout = row(row(), name=title_prefix + ' amino acid pauses')

    for key in amino_acid_df_dict.keys():
        print(key)
        key_title = get_key_title(title_prefix, key)
        amino_acid_df = amino_acid_df_dict.get(key)
        p = figure(title=key_title,
                   x_axis_label="distance from amino acid", y_axis_label="5'seq read counts")

        for aa in amino_acid_df.index:
            c = codons.Codons.AMINO_ACID_COLORS.get(aa)
            line = p.line(amino_acid_df.columns, amino_acid_df.loc[aa,], legend=aa,
                          line_color=c)
            hover = HoverTool(tooltips=[('distance', '@x'), ('count', '@y')], renderers=[line])
            p.add_tools(hover)

        p.legend.click_policy = "hide"
        mainLayout.children[0].children.append(p)
        export_images(p, key_title, png_dir, svg_dir)

    return mainLayout


def bokeh_aa_pause_linechart(title, amino_acid_df):
    p = figure(title=title,
               x_axis_label="distance from amino acid", y_axis_label="5'seq read counts")

    for aa in amino_acid_df.index:
        c = codons.Codons.AMINO_ACID_COLORS.get(aa)
        # line = p.line(amino_acid_df.columns, amino_acid_df.loc(aa, ), legend=aa, line_color=RGB(c[0], c[1], c[2]))
        line = p.line(amino_acid_df.columns, amino_acid_df.loc(aa, ), legend=aa, line_color=c)
        hover = HoverTool(tooltips=[('distance', '@x'), ('count', '@y')], renderers=[line])
        p.add_tools(hover)

    p.legend.click_policy = "hide"
    return p


def bokeh_frame_line_charts(title_prefix, region_count_series_dict, w=20,
                            lib_size_dict=None, png_dir=None, svg_dir=None):
    if region_count_series_dict is None:
        return None

    #logging.getLogger(config.FIVEPSEQ_LOGGER).info(
    #    "Making frame line charts %s with options: smoothing window(%d), " %
    #    (title_prefix, w))

    my_y_label = "5' counts per base (averaged over %d codon window), RPM" % w

    figure_list = []

    frame_cols = {"f0": "#ee6677", "f1": "#5ab76f", "f2": "#66ccee"}

    key0 = list(region_count_series_dict)[0]

    for sample_ind in range(len(region_count_series_dict.get(key0))):
        sample = list(region_count_series_dict.get(key0))[sample_ind]
        sample_title = get_key_title(title_prefix, sample)
        sample_charts = row(row(), name=sample_title)

        for region in [FivePSeqCounts.START, FivePSeqCounts.MID, FivePSeqCounts.TERM]:
            my_x_label = "distance from %s" % region

            p = figure(title=sample_title + "_" + region, x_axis_label=my_x_label, y_axis_label=my_y_label)

            legend_items = []

            count_series = region_count_series_dict.get(region).get(sample)

            if count_series is None:
                # return None
                continue

            try:
                for f in (0, 1, 2):
                    c = frame_cols.get("f" + str(f))
                    # adjust the frame so that the position 0 corresponds to frame 0, 1 to 1 and 2 to 2.
                    f_start_ind = list(np.where([d % 3 == f for d in count_series.D][0:3])[0])[0]
                    f_counts = count_series.C[f_start_ind::3]
                    f_rpm_counts = [10 ** 6 * x / lib_size_dict.get(sample) for x in f_counts]
                    y = average_counts_over_window(f_rpm_counts, w)
                    source = ColumnDataSource(
                        data=dict(
                            x=list(count_series.D[f_start_ind::3]),
                            y=list(y),
                            name=["f" + str(f)] * len(y)
                        ))
                    line = p.line('x', 'y', line_color=c, line_width=2, source=source)
                    legend_items.append((sample + ": F" + str(f), [line]))
                    hover = HoverTool(tooltips=[('name', '@name'), ('position', '@x')], renderers=[line])
                    p.add_tools(hover)

            except Exception as e:
                print("Error at key %s, reason; %s" % (sample, str(e)))
                return None

            legend = Legend(items=legend_items, location=LEGEND_ALIGNMENT)
            p.add_layout(legend, LEGEND_POSITION)
            p.legend.click_policy = "hide"

            sample_charts.children[0].children.append(p)

            # figures for exporting
            export_images(p, sample_title + "-" + region, png_dir, svg_dir)

        figure_list.append(sample_charts)

    frame_grid = gridplot(figure_list, ncols=1)

    return frame_grid


def average_counts_over_window(counts, w):
    """
    Counts per position averaged over the given window size. The first positions are not averaged.

    :param counts: vector of counts
    :param w: window size
    :return: vector of average counts (the same length as the original vector)
    """

    av_vec = [0] * len(counts)
    for i in range(len(counts)):
        if i < w / 2:
            if len(counts) > 2 * i:
                av_vec[i] = sum(counts[:2 * i + 1]) / (2 * i + 1)
            else:
                av_vec[i] = sum(counts) / len(counts)
        else:
            if len(counts) >= i + int(np.ceil(w / 2)):
                av_vec[i] = sum(counts[i - int(np.floor(w / 2)): i + int(np.ceil(w / 2))]) / w
            else:
                av_vec[i] = sum(counts[i:]) / (len(counts) - i)

    return av_vec


"""
def bokeh_pca_plot(title, codon_df_dict, color_dict, start_pos = '-20', end_pos = '-2',
                   png_dir=None, svg_dir=None):

    def scale(df):
        sdf = df.copy(deep = True)
        for i in range(sdf.shape[0]):
            sdf.iloc[i, :] /= sum(sdf.iloc[i, :]) + 1
        return sdf

    def mat_to_vec(sdf):
        vec = []
        vecnames = []
        for j in range(sdf.shape[0]):
            start_ind = [i for i, j in enumerate(sdf.columns) if j == start_pos][0]
            stop_ind = [i for i, j in enumerate(sdf.columns) if j == end_pos][0]
            vec += list(sdf.iloc[j, start_ind:stop_ind])
            vecnames += [sdf.index[j] + codon for codon in list(sdf.columns[start_ind:stop_ind])]
        return vec, vecnames

    vec_list = []
    vecnames = None
    colors = []
    names = []
    for s in codon_df_dict.keys():
        colors.append(color_dict.get(s))
        names.append(s)
        vec, vecnames = mat_to_vec(scale(codon_df_dict[s]))
        vec_list.append(vec)

    data = pd.DataFrame(vec_list, columns=vecnames)
    pca = PCA(n_components=2)

    principalComponents = pca.fit_transform(data)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

    source = ColumnDataSource(data = dict(
        x = principalDf['PC1'],
        y = principalDf['PC2'],
        col = colors,
        sample = names
    ))

    # draw
    pc1_var = list(pca.explained_variance_ratio_)[0] * 100
    pc2_var = list(pca.explained_variance_ratio_)[1] * 100

    p = figure(title=title, height=500,
               x_axis_label="PC1 (%.2f %% variance)" % pc1_var,
               y_axis_label="PC2 (%.2f %% variance)" % pc2_var)

    circles = p.circle('x', 'y', color='col', legend='sample',
                       fill_alpha=0.8, size=20, line_color=None,
                       source=source)
    hover = HoverTool(tooltips = [('sample', '@sample')], renderers = [circles])
    p.add_tools(hover)
    p.legend.click_policy = "hide"

    if png_dir is not None:
        p_png = figure(title = title, height = 500)
        p_png.circle('x', 'y', color=colors, legend_group='sample',
                       fill_alpha=0.8, size=20, line_color=None,
                       source=source)
        export_images(p_png, title, png_dir=png_dir)

    if svg_dir is not None:
        p_svg = figure(title = title, height = 500)
        p_svg.circle('x', 'y', color=colors, legend_group='sample',
                       fill_alpha=0.8, size=20, line_color=None,
                       source=source)
        export_images(p_svg, title, svg_dir=svg_dir)
    return p
"""


def export_images(p, title, png_dir=None, svg_dir=None):
    if png_dir is not None:
        try:
            png_f = os.path.join(png_dir, title + ".png")
            p.background_fill_color = None
            p.border_fill_color = None
            export_png(p, filename=png_f)
        except Exception as e:
            logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Problem exporting figure %s. Reason: %s"
                                                              % (title, str(e)))
    if svg_dir is not None:
        try:
            svg_f = os.path.join(svg_dir, title + ".svg")
            p.output_backend = "svg"
            export_svgs(p, filename=svg_f)
        except Exception as e:
            logging.getLogger(config.FIVEPSEQ_LOGGER).warning("Problem exporting figure %s. Reason: %s"
                                                              % (title, str(e)))


def get_key_title(title, key, scale=False):
    return key + "-" + title


def get_random_color():
    return cl.scales['9']['qual']['Set3'][random.randint(0, 8)]