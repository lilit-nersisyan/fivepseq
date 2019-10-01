Welcome to fivepseq readme!
--------
This program reports 5'-footprint properties from 5Pseq reads in a given alignment file. 
 
There are to modes to run fivepseq: 
count - generates count vectors and meta-vectors from the given alignment file and stores them in the output directory as text files
plot - generates plots based on input count data folder and stores those to an html file
count_and_plot - performs counts and generates plots afterwards

Preprocessing from FASTQ files (!NEW)
--------
A script for processing raw fastq files from yeast has been added to preprocess_scripts directory.
It takes raw fastq files and

- does quality checks (with FASTQC and MULTIQC),

- trims adapters,

- extracts UMI,

- generates a STAR index if not provided

- maps with STAR

- deduplicates reads with umitools

- analyzes RNA transcript content (relative content of coding versus non-coding RNA). The final stats are in the align_rna/rna_stats.txt


*Run the script with options:*


-f directory with fastq files (leave only the files you'd like to use (Read1 for fivepseq))

-o output directory

-g genome (fasta) file path

-a annotation (gff) file path

-i STAR index path, if you'd like to use existing index

-s specify which steps of the pipeline you'd like to skip. Possible values are:

   c   skip trimming adapters with cutadapt

   u   skip UMI extraction

   d   skip deduplication after alignment

   q   skip quality check: FASTQC and MULTIQC

   m   skip mapping

   or any combination of these characters, e.g. use -s cudqm to skip all

Before you start fivepseq
--------
- Make sure you're using python version 2.7
- Install the latest stable version of fivepseq by
- - cd /proj/sllstore2017018/lilit/fivepseq_latest:
- - python setup.py install

Usage
--------

fivepseq \\

-g path_to_your_genome_fasta \\

-a path_to_your_annotation_gff (sorry no gtf for now) \\

-b path_to_one_or_many_bam_files (many files should be provided with pattern (e.g. "parent_dir/*.bam") within brackets) \\

-o output_directory \\

-t title_of_the_run \\

[optional arguments] \\



optional argument: --conflicts

The conflict mode specifies how to deal with files/folders that alraedy exist. There are three options you may choose from:

- add (default) - only adds missing files to existing directory

- overwrite - overwrites all the files in the existing (in case) output directory

- alt_dir - uses alternative directory by appending '+' suffix to existing (in case) output directory


optional argument: --op (default = 0)

This arguments sets the p value threshold for outlier detection: point with less than the --op probability of
falling into Poisson distribution will be down-sampled. If you want to turn off downsampling, set the --op to -1.


optional argument: --ds

A constant value for down-sampling. Instead of outlier detection, values less than this constant will be down-sampled
to match --ds.


optional argument: -gf/-genefilter

Supply a text file with newline-separated list of gene ids you'd like to filter/use. The names should correspond to those present under the gene_id tag in the gff file.
Note, only these genes will be used in all the calculations.


optional argument: -gs/-geneset
This option provides a possibility to compared plots for different samples. Supply a tab separated text file, with the following structure:
Column names: gene_attribute (e.g. Name)->geneset
Rows: value_of_the_attribute->geneset_name

Note, the gene_attribute is the attribute name in the gtf or gff file. In case of gff, the attribute in the cds feature will be considered.
With this option, fivepseq will generate a separate plotting directory called genesets, with tabbed-plots to either compare samples for each geneset, or genesets for each sample.
The counts folder will also be divided according to the geneset used. The default folder will be named protein_coding.


optional argument: --loci-file

This option requires a file with coordinates of the loci (e.g. RBP binding coordinates), relative to which, the user wants to generate scatter-plots.
The file should be tab-separated, with the following structure:
Columns: chr->str->start->end->symbol
Rows: the chromosome name, strand (+ or -), start and end coordinates and the name of the RBP (or the locus).

This feature is in beta-, so scatter plots of reads relative to all the loci combined will be plotted underneath the main html file.
Four different plots will correspond to reads located in (1) 3UTR and CDS regions, (2) only 3UTR, (3) only 5UTR and (4) only CDS.


Note!
-------
bai index files should be in the same directory as the bam files


Have fun!
