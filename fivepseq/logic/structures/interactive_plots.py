import plotly
import plotly.plotly as py
import plotly.graph_objs as go

import pandas as pd

plotly.tools.set_credentials_file(username='lilit_nersisyan', api_key='OSKYQ02w2w3f34hJA55M')
# df = pd.read_csv("/proj/sllstore2017018/lilit/5pseq_microbiome/fivepseq_out_l_pla/meta_counts_TERM.txt", sep = "\t")
# df = pd.read_csv("/proj/sllstore2017018/lilit/5pseq_human/fivepseq_out1/meta_counts_TERM.txt", sep = "\t")
#df = pd.read_csv("/proj/sllstore2017018/lilit/5pseq_human/fivepseq_UMI_Hela-rep1_S7_R1_001/meta_counts_TERM.txt", sep = "\t")
df = pd.read_csv("/proj/sllstore2017018/lilit/5pseq_human/fivepseq_HelaFrag-rep1/meta_counts_TERM.txt", sep = "\t")


trace_high = go.Scatter(
    x=df.D1,
    y=df['Count'],
    name = "Count",
    line = dict(color = '#17BECF'),
    opacity = 0.8)

# trace_low = go.Scatter(
#     x=df.Date,
#     y=df['AAPL.Low'],
#     name = "AAPL Low",
#     line = dict(color = '#7F7F7F'),
#     opacity = 0.8)

data = [trace_high]

layout = dict(
    title='Meta counts on terminus +/- 100 bp',
    xaxis=dict(
        rangeselector=dict(
            buttons=list([
                dict(count=1,
                     label='1bp'),
                dict(step='all')
            ])
        ),
        rangeslider=dict(
            visible = True
        )
    )
)

fig = dict(data=data, layout=layout)
py.iplot(fig, filename = "Meta counts on spanning terminus HelaFrag-rep1")