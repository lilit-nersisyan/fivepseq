import os

import pandas as pd
import numpy as np


dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"
suffix = "counts_FULL_LENGTH.txt"
hela_rep1_full = pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", suffix),
        sep="\t")
hela_chx_rep1_full = pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaCHX-rep1", suffix),
        sep="\t")
hela_frag_rep1_full = pd.read_csv(
        os.path.join(dir_5pseq_human, "fivepseq_HelaFrag-rep1", suffix),
        sep="\t")


diff_untr = hela_rep1_full - hela_frag_rep1_full
diff_chx = hela_chx_rep1_full - hela_frag_rep1_full

ind = diff_untr > 0 and diff_chx > diff_untr