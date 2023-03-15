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
    START_CODON = "ATG"

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

    AMINO_ACID_SINGLELETTER_TABLE = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLY": "G",
        "GLU": "E",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "TERM": "*"
    }

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

    @staticmethod
    def get_dicodon_table():
        dicodon_table = {}
        for codon1 in Codons.CODON_TABLE.keys():
            for codon2 in Codons.CODON_TABLE.keys():
                dicodon_table.update(
                    {codon1 + codon2: Codons.CODON_TABLE.get(codon1) + '-' + Codons.CODON_TABLE.get(codon2)})
        return dicodon_table

    @staticmethod
    def get_tricodon_table():
        tricodon_table = {}
        for codonE in Codons.CODON_TABLE.keys():
            for codonP in Codons.CODON_TABLE.keys():
                for codonA in Codons.CODON_TABLE.keys():
                    tricodon_table.update({codonE + codonP + codonA:
                                               Codons.get_peptide_from_codon_list([codonE, codonP, codonA])
                                           })
        return tricodon_table

    @staticmethod
    def get_peptide_from_codon_list(codon_list):
        return '-'.join([Codons.CODON_TABLE[codon] for codon in codon_list]) +\
               ' (' + ''.join([Codons.AMINO_ACID_SINGLELETTER_TABLE[Codons.CODON_TABLE[codon]] for codon in codon_list]) +\
               ')'


    @staticmethod
    def get_dipeptide_list():
        dipeptide_list = list()
        for aaP in Codons.AMINO_ACID_TABLE.keys():
            for aaA in Codons.AMINO_ACID_TABLE.keys():
                dipeptide_list.append(
                    '-'.join([aaP, aaA]) +
                    ' (' +
                    ''.join([Codons.AMINO_ACID_SINGLELETTER_TABLE.get(aaP),
                             Codons.AMINO_ACID_SINGLELETTER_TABLE.get(aaA)
                             ]) +
                    ')'
                )

        return dipeptide_list

    @staticmethod
    def get_tripeptide_list():
        tripeptide_list = list()
        for aaE in Codons.AMINO_ACID_TABLE.keys():
            for aaP in Codons.AMINO_ACID_TABLE.keys():
                for aaA in Codons.AMINO_ACID_TABLE.keys():
                    tripeptide_list.append(
                        '-'.join([aaE, aaP, aaA]) +
                        ' (' +
                        ''.join([Codons.AMINO_ACID_SINGLELETTER_TABLE.get(aaE),
                                 Codons.AMINO_ACID_SINGLELETTER_TABLE.get(aaP),
                                 Codons.AMINO_ACID_SINGLELETTER_TABLE.get(aaA)
                                 ]) +
                        ')'
                    )

        return tripeptide_list

    @staticmethod
    def get_dicodon_full_index():
        dicodon_index = list()
        for codonP in Codons.CODON_TABLE.keys():
            for codonA in Codons.CODON_TABLE.keys():
                dicodon_index.append(
                    '-'.join([Codons.CODON_TABLE.get(codonP),
                              Codons.CODON_TABLE.get(codonA)]) +
                    '_' +
                    ''.join(codonP + codonA))
        return dicodon_index

    @staticmethod
    def get_tricodon_full_index():
        tricodon_index = list()
        for codonE in Codons.CODON_TABLE.keys():
            for codonP in Codons.CODON_TABLE.keys():
                for codonA in Codons.CODON_TABLE.keys():
                    tricodon_index.append(
                        '-'.join([Codons.CODON_TABLE.get(codonE),
                                  Codons.CODON_TABLE.get(codonP),
                                  Codons.CODON_TABLE.get(codonA)]) +
                        '_' +
                        ''.join(codonE + codonP + codonA))
        return tricodon_index

    def get_color20(i):
        tab20 = cm.get_cmap("tab20", 20)
        col = tab20(range(20))[i] * 255
        rgb = RGB(col[0], col[1], col[2])
        return rgb

    AMINO_ACID_COLORS = dict(zip(AMINO_ACID_TABLE.keys(), map(get_color20, list(range(20)))))
