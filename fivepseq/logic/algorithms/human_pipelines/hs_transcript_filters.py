import os

import pandas as pd

dir_5pseq_human = "/proj/sllstore2017018/lilit/5pseq_human"

transcript_assembly = pd.read_csv(
    os.path.join(dir_5pseq_human, "fivepseq_Hela-rep1", "transcript_assembly.txt"),
    sep="\t")



