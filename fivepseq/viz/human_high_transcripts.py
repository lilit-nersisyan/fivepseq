import os

import pandas as pd
import colorlover as cl
import numpy as np
from fivepseq.viz.bokeh_plots import bokeh_transcript_scatter_plot

from fivepseq.logic.structures.fivepseq_counts import CountManager, FivePSeqCounts

dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"

transcript_assembly = pd.read_csv(
    os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "transcript_assembly.txt"),
    sep="\t")

hela_chx_counts = CountManager.read_counts_as_list(
    os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", "counts_FULL_LENGTH.txt"))


high_transcripts = pd.read_csv("/proj/sllstore2017018/lilit/5pseq_human/resources/top1000.transcripts.txt")
t_ids = transcript_assembly.iloc[:,0]
t_ids = [w.replace("transcript:", "") for w in t_ids]
chx_high_ind = [i for i, item in enumerate(t_ids) if item in high_transcripts]

bokeh_transcript_scatter_plot("HelaCHX-rep1",
                              {"HelaCHX-rep1": hela_chx_counts},
                              transcript_assembly,
                              {"HelaCHX-rep1" : cl.to_numeric(cl.scales['9']['qual']['Set3'])[2]},
                              FivePSeqCounts.TERM, 500,
                              index_filter=chx_high_ind, min_count=0)
