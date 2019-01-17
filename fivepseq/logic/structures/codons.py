"""
This class keeps codon annotations as public variables to be referenced from other functions.
"""
from bokeh.colors import RGB
from matplotlib import cm


class Codons:
    def __init__(self):
        pass

    STOP_AMBER = "TAG"
    STOP_OCHRE = "TAA"
    STOP_UMBER = "TGA"
    stop_codons = [STOP_AMBER, STOP_OCHRE, STOP_UMBER]

    AMINO_ACID_TABLE = {
        "ALA": ["GCT", "GCC", "GCA", "GCG"],
        "ARG": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "ASN": ["AAT", "AAC"],
        "ASP": ["GAT", "GAC"],
        "CYS": ["TGT", "TGC"],
        "GLN": ["CAA", "CAG"],
        "GLU": ["GAA", "GAG"],
        "GLY": ["GGT", "GGC", "GGA", "GGG"],
        "HIS": ["CAT", "CAC"],
        "ILE": ["ATT", "ATC", "ATA"],
        "LEU": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "LYS": ["AAA", "AAG"],
        "PHE": ["TTT", "TTC"],
        "PRO": ["CCT", "CCC", "CCA", "CCG"],
        "SER": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "THR": ["ACT", "ACC", "ACA", "ACG"],
        "TRP": ["TGG"],
        "TYR": ["TAT", "TAC"],
        "VAL": ["GTT", "GTC", "GTA", "GTG"],
        "TERM": ["TAA", "TGA", "TAG"],
        "MET": ["ATG"]}

    CODON_TABLE = {
        "GCT": "ALA",
        "GCC": "ALA",
        "GCA": "ALA",
        "GCG": "ALA",
        "CGT": "ARG",
        "CGG": "ARG",
        "CGC": "ARG",
        "CGA": "ARG",
        "AGA": "ARG",
        "AGG": "ARG",
        "AAT": "ASN",
        "AAC": "ASN",
        "GAT": "ASP",
        "GAC": "ASP",
        "TGT": "CYS",
        "TGC": "CYS",
        "CAA": "GLN",
        "CAG": "GLN",
        "GAA": "GLU",
        "GAG": "GLU",
        "GGT": "GLY",
        "GGC": "GLY",
        "GGA": "GLY",
        "GGG": "GLY",
        "CAT": "HIS",
        "CAC": "HIS",
        "ATT": "ILE",
        "ATC": "ILE",
        "ATA": "ILE",
        "TTA": "LEU",
        "TTG": "LEU",
        "CTT": "LEU",
        "CTC": "LEU",
        "CTA": "LEU",
        "CTG": "LEU",
        "AAA": "LYS",
        "AAG": "LYS",
        "TTT": "PHE",
        "TTC": "PHE",
        "CCT": "PRO",
        "CCC": "PRO",
        "CCA": "PRO",
        "CCG": "PRO",
        "TCT": "SER",
        "TCC": "SER",
        "TCA": "SER",
        "TCG": "SER",
        "AGT": "SER",
        "AGC": "SER",
        "ACT": "THR",
        "ACC": "THR",
        "ACA": "THR",
        "ACG": "THR",
        "TGG": "TRP",
        "TAT": "TYR",
        "TAC": "TYR",
        "GTT": "VAL",
        "GTC": "VAL",
        "GTA": "VAL",
        "GTG": "VAL",
        "TAA": "TERM",
        "TGA": "TERM",
        "TAG": "TERM",
        "ATG": "MET"}

    def get_color20(i):
        tab20 = cm.get_cmap("tab20", 20)
        col = tab20(range(20))[i] * 255
        rgb = RGB(col[0], col[1], col[2])
        return rgb

    AMINO_ACID_COLORS = dict(zip(AMINO_ACID_TABLE.keys(), map(get_color20, list(range(20)))))
