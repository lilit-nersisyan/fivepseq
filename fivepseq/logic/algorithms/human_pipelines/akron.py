import time

from fivepseq.logic.structures.counts import FivePSeqCounts

from fivepseq.util.readers import BamReader, AnnotationReader, FastaReader

from fivepseq import config

from fivepseq.__main__ import FivepseqArguments, setup_output_dir, setup_logger

#fivepseq_arguments = FivepseqArguments()
#print "Call to fivepseq with arguments:"
#fivepseq_arguments.print_args()

# startup
#setup_output_dir()
#setup_logger()
#config.logger.info("Fivepseq started")
#start_time = time.clock()

bam = "/proj/sllstore2017018/lilit/akron/riboseq/staralign.accurate2/SRR3306589_1.fastqAligned.sortedByCoord.out.bam"
annot = "/proj/sllstore2017018/lilit/5pseq_human/gff/Homo_sapiens.GRCh38.93.gff3"
genome = "/proj/sllstore2017018/lilit/5pseq_human/genome2/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"

bam_reader = BamReader(bam)
annotation_reader = AnnotationReader(annot, 20000)
fasta_reader = FastaReader(genome)
fivepseq_counts = FivePSeqCounts(bam_reader.alignment, annotation_reader.annotation, fasta_reader.genome,
                                 0)


##############
# given the coordinates of riboseq summits, find the distances from these coordinates.
# Use the count files generated from the whole genome.


