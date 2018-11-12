import os

import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import colorlover as cl
import pandas as pd
import numpy as np

from fivepseq.viz.plotly_plots import PlotlyPlots


class MetaCountPlots(PlotlyPlots):

    def __init__(self):
        PlotlyPlots.__init__(self)

    def plot_meta_counts(self, meta_count_dict, color_dict, filename):
        """

        :param self:
        :param meta_count_dict: {sample:[meta_counts_df]}
        :param color_dict: {sample:color}
        :param filename:
        :return:
        """
        self.__init__()
        for i in range(0, len(meta_count_dict) - 1):
            name1 = meta_count_dict.keys()[i]
            name2 = meta_count_dict.keys()[i + 1]
            if not all(meta_count_dict.get(name1).D1 == meta_count_dict.get(name2).D1):
                raise ValueError("indices of df %d and df %d did not match" % (i, i + 1))

        data = []
        for name in meta_count_dict.keys():
            meta_df = meta_count_dict.get(name)
            trace = go.Scatter(
                x=meta_df.D1,
                y=meta_df.Count,
                name=name,
                line=dict(color=color_dict.get(name)),
                opacity=0.8
            )
            data.append(trace)

        layout = dict(
            title="Meta counts terminal region",
            xaxis=dict(
                rangeselector=dict(
                    buttons=list([
                        dict(step="all")
                    ])
                ),
                rangeslider=dict(
                    visible=True
                )
            )
        )

        fig = dict(data=data, layout=layout)
        py.iplot(fig, filename=filename)

