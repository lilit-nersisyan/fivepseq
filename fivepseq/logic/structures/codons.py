"""
This class keeps codon annotations as public variables to be referenced from other functions.
"""


class Codons:
    def __init__(self):
        pass

    STOP_AMBER = "TAG"
    STOP_OCHRE = "TAA"
    STOP_UMBER = "TGA"
    stop_codons = [STOP_AMBER, STOP_OCHRE, STOP_UMBER]
