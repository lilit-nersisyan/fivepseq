from Bio import SeqIO


class Genome:
    fasta_file = None
    genome_dict = None

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.genome_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file), "fasta"))
