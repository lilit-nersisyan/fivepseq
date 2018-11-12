import fivepseq
import os

import plotly
import plotly.plotly as py
import plotly.io as pio
import plotly.graph_objs as go
import colorlover as cl
import pandas as pd
import numpy as np
from fivepseq.logic.structures import fivepseq_counts

from fivepseq import config


class InteractivePlots:

    def __init__(self):
        plotly.tools.set_credentials_file(username='lilit_nersisyan', api_key='OSKYQ02w2w3f34hJA55M')

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

    def triangle_plots(self, frame_counts_dict, color_dict, filename, max_threshold, min_threshold):

        raw_data = []
        max_count = 0
        for specie in frame_counts_dict.keys():
            df = frame_counts_dict.get(specie)
            if max_threshold is not None:
                # filter rows with counts more than max_threshold
                df = df[(df[['F0', 'F1', 'F2']] <= max_threshold).all(axis=1) & (df[['F0', 'F1', 'F2']] >= min_threshold).all(axis=1)]
                print "Transcript count for %s is %d, filtered count: %d" % (specie, len(df), len(df))

            df['color'] = color_dict.get(specie)
            raw_data = raw_data + df.to_dict(orient='index').values()
            max_value = max(df[['F0','F1','F2']].max())
            if max_value > max_count:
                max_count = max_value
        """
        rawData = [
            {'journalist': 75, 'developer': 25, 'designer': 0, 'label': 'point 1'},
            {'journalist': 70, 'developer': 10, 'designer': 20, 'label': 'point 2'},
            {'journalist': 75, 'developer': 20, 'designer': 5, 'label': 'point 3'},
            {'journalist': 5, 'developer': 60, 'designer': 35, 'label': 'point 4'},
            {'journalist': 10, 'developer': 80, 'designer': 10, 'label': 'point 5'},
            {'journalist': 10, 'developer': 90, 'designer': 0, 'label': 'point 6'},
            {'journalist': 20, 'developer': 70, 'designer': 10, 'label': 'point 7'},
            {'journalist': 10, 'developer': 20, 'designer': 70, 'label': 'point 8'},
            {'journalist': 15, 'developer': 5, 'designer': 80, 'label': 'point 9'},
            {'journalist': 10, 'developer': 10, 'designer': 80, 'label': 'point 10'},
            {'journalist': 20, 'developer': 10, 'designer': 70, 'label': 'point 11'},
        ]
        """

        def makeAxis(title, tickangle):
            return {
                'title': title,
                'titlefont': {'size': 20},
                'tickangle': tickangle,
                'tickfont': {'size': 15},
                'tickcolor': 'rgba(0,0,0,0)',
                'ticklen': 5,
                'showline': True,
                'showgrid': True
            }

        data = [{
            'type': 'scatterternary',
            'mode': 'markers',
            'a': [i for i in map(lambda x: x['F0'], raw_data)],
            'b': [i for i in map(lambda x: x['F1'], raw_data)],
            'c': [i for i in map(lambda x: x['F2'], raw_data)],
            'marker': {
                'symbol': 100,
                'size': 14,
                'color': [i for i in map(lambda x: x['color'], raw_data)],
                'line': {'width': 2}
            },
        }]

        layout = {
            'ternary': {
                'sum': max_count,
                'aaxis': makeAxis('F0', 0),
                'baxis': makeAxis('F1', 45),
                'caxis': makeAxis('F2', -45)
            },
            'annotations': [{
                'showarrow': False,
                'text': filename,
                'x': 0.5,
                'y': 1.3,
                'font': {'size': 15}
            }]
        }

        fig = {'data': data, 'layout': layout}
       # FIXME problem with orca, orca executable not found, will not save image: pio.write_image(fig,filename + ".jpg")
        py.iplot(fig, validate=False, filename=filename)

    def transcript_filter(self, df_dict, f0_filter, f1_filter, f2_filter):
        t_inds = []
        for key in df_dict.keys():
            df = df_dict.get(key)
            t_ind = [True] * len(df)
            if f0_filter is not None:
                if f0_filter[0] == 'max':
                    t_ind = np.logical_and(t_ind, df['F0'] < f0_filter[1])
                else:
                    t_ind = np.logical_and(t_ind, df['F0'] > f0_filter[1])
            if f1_filter is not None:
                if f1_filter[0] == 'max':
                    t_ind = np.logical_and(t_ind, df['F1'] < f1_filter[1])
                else:
                    t_ind = np.logical_and(t_ind, df['F1'] > f1_filter[1])
            if f2_filter is not None:
                if f2_filter[0] == 'max':
                    t_ind = np.logical_and(t_ind, df['F2'] < f2_filter[1])
                else:
                    t_ind = np.logical_and(t_ind, df['F2'] > f2_filter[1])
            t_inds.append(t_ind)
        t_index_dict = dict(zip(df_dict.keys(), t_inds))
        return t_index_dict

    def filtered_meta_counts(self, t_ind, dir, fivepseq_prefix="fivepseq_", count_suffix="count_TERM.txt"):
        meta_count_dict = {}
        for sample in t_ind.keys():
            count_file = os.path.join(dir, fivepseq_prefix + sample, count_suffix)
            count_df = pd.read_csv(count_file, header=None, sep='\t')
            count_df = count_df[t_ind.get(sample)]
            counts = count_df.values
            #FIXME temporarily added 1000
            meta_counts = fivepseq.logic.structures.fivepseq_counts.FivePSeqCounts.compute_meta_counts(counts, 1000)
            meta_count_dict[sample] = meta_counts
        return meta_count_dict

##################################################################################
#                               TEST PLOTS: SCATTERPLOTS                         #
##################################################################################


dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"
meta_count_rep1_term_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HEK293-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HEK293CHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293CHX-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HEK293Frag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293Frag-rep1", "meta_counts_TERM.txt"),
        sep="\t")
}


meta_count_Hela_rep1_term_dict = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_TERM.txt"),
        sep="\t"),
}

meta_count_dict_rep2 = {
    "Hela-rep2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep2", "meta_counts_TERM.txt"),
        sep="\t"),
    "HelaCHX-rep2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep2", "meta_counts_TERM.txt"),
        sep="\t"),
    "HelaFrag-rep2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep2", "meta_counts_TERM.txt"),
        sep="\t"),
    "HEK293-rep2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293-rep2", "meta_counts_TERM.txt"),
        sep="\t"),
    "HEK293CHX-rep2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293CHX-rep2", "meta_counts_TERM.txt"),
        sep="\t"),
    "HEK293Frag-rep2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293Frag-rep2", "meta_counts_TERM.txt"),
        sep="\t"),
    "SRR7902417_1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_SRR7902417_1", "meta_counts_TERM.txt"),
        sep="\t"),
    "SRR7902417_2": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_SRR7902417_2", "meta_counts_TERM.txt"),
        sep="\t")
}
colors_dict = dict(zip(meta_count_rep1_term_dict.keys(), cl.scales['6']['qual']['Set1']))
colors_dict_rep2 = dict(zip(meta_count_dict_rep2.keys(), cl.scales['8']['qual']['Dark2']))

meta_count_dict_start_rep1 = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "meta_counts_START.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "meta_counts_START.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "meta_counts_START.txt"),
        sep="\t"),
    "HEK293-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293-rep1", "meta_counts_START.txt"),
        sep="\t"),
    "HEK293CHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293CHX-rep1", "meta_counts_START.txt"),
        sep="\t"),
    "HEK293Frag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HEK293Frag-rep1", "meta_counts_START.txt"),
        sep="\t")
}
colors_dict_start_rep1 = dict(zip(meta_count_dict_start_rep1.keys(), cl.scales['6']['qual']['Dark2']))

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
print "here"

##################################################################################
#                               TEST PLOTS: TRIANGLEPLOTS                        #
##################################################################################

dir_5pseq_microbiome = "/proj/sllstore2017018/lilit/5pseq_microbiome"
frame_counts_dict_lpla = {
    "lpla_untreated": pd.read_csv(
        os.path.join(dir_5pseq_microbiome, "fivepseq_lpla_untreated", "frame_counts_START.txt"),
        sep="\t"),
    "lpla_fragmented": pd.read_csv(
        os.path.join(dir_5pseq_microbiome, "fivepseq_l_pla_fragmented", "frame_counts_TERM.txt"),
        sep="\t")
}

frame_count_dict_start_rep1 = {
    "Hela-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "frame_counts_START.txt"),
        sep="\t"),
    "HelaCHX-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "frame_counts_START.txt"),
        sep="\t"),
    "HelaFrag-rep1": pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", "frame_counts_START.txt"),
        sep="\t"),
}

frame_count_dict_term_rep1 = {
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

interactive_plots = InteractivePlots()
#t_ind = interactive_plots.transcript_filter(frame_count_dict_term_rep1, f0_filter=['max', 20], f1_filter=['max', 40], f2_filter=['min', 60])
#filtered_meta_count_dict = interactive_plots.filtered_meta_counts(t_ind, dir_5pseq_human, count_suffix="count_TERM.txt")
#interactive_plots.plot_meta_counts(meta_count_dict=filtered_meta_count_dict, color_dict=colors_dict, filename="Terminal counts, human, rep1, f2_filters")
# interactive_plots.plot_meta_counts(meta_count_dict_rep2, colors_dict_rep2, "Terminal counts, human, rep2")
# interactive_plots.plot_meta_counts(meta_count_dict_start_rep1, colors_dict_start_rep1, "Start counts, human, rep1")
# interactive_plots.plot_meta_counts(meta_count_dict_lpla, colors_dict_lpla, "Terminal counts, L_pla")
#interactive_plots.triangle_plots(frame_counts_dict_lpla, colors_dict_lpla, "Triangle plots for l_pla")
#interactive_plots.triangle_plots(frame_count_dict_start_rep1, colors_dict_start_rep1, "Triangle plots for Hela rep-1")
interactive_plots.triangle_plots(frame_count_dict_term_rep1, colors_dict_start_rep1, "Triangle plots for Hela rep-1, TERM", max_threshold=500, min_threshold = 30)
