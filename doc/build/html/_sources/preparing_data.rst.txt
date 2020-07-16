.. _preparing_data:



Preparing the data
==================

Required files
-----------------

In order to run fivepseq, the following files are required:

|    Aligned reads (**.bam**)
|    Alignment index (**.bai**)
|    Genomic sequence file (**.fasta** / **.fa**)
|    Genomic annotation file (**.gff/ .gtf**)

This section is to help you generate alignment bam files from your fastq files. If you already have alignment, genome and annotation files you may skip this section and proceed to fivepseq downstream analysis. 


Finding genome and annotation files
--------------------------------------

There are multiple resources which maintain and distribute genome sequence and annotation data, for instance `Ensembl <https://www.ensembl.org/index.html>`_, `RefSeq <https://www.ncbi.nlm.nih.gov/assembly/>`_, and `SGD <https://www.yeastgenome.org/>`_.
However, the naming conventions and format may differ slightly depending on the resource. 
For this reason, your genomic sequence and annotation files should be retrieved from the same database.

**!!NOTE!!** If you start from alignment files, make sure that the annotation and genome files you provide to fivepseq for downstream analysis are the same as those used during alignment. 


Preprocessing fastq files
-----------------------------------

Fastq files need to be preprocessed and aligned to the reference genome before proceeding to fivepseq downstream analysis. Preprocessing proceeds with the following steps: 

- quality checks (with FASTQC and MULTIQC),
- adapter and quality based trimming,
- UMI extraction (if the library was generated with UMIs),
- mapping to reference
- read deduplication (if the library was generated with UMIs),
- generation of bedgraph files to visualize 5'P counts in genome viewer

An example of pre\-processing pipeline can be found in the **preprocess_scripts/fivepseq_preprocess.sh** directory in your fivepseq installation folder, and can also be downloaded from `here <https://github.com/lilit-nersisyan/fivepseq/blob/master/preprocess_scripts/fivepseq_preprocess.sh>`_. 

In order to run this pipeline, you need to have access to common bioinformatics software such as `STAR <https://github.com/alexdobin/STAR>`_, `UMI-tools <https://github.com/CGATOxford/UMI-tools>`_, `bedtools <https://bedtools.readthedocs.io/en/latest/>`_, `Samtools <http://www.htslib.org/>`_, `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_, `MultiQC <https://multiqc.info/>`_ and `cutadapt <https://github.com/marcelm/cutadapt>`_. 

To use it, navigate to the directory where the script is located and use the following command in the prompt: 

.. code-block:: shell

 ./fivepseq_preprocess.sh \
   -f [path to directory containing fastq files] \
   -g [path to genome fasta] \
   -a [path to annotation gff/gtf] \
   -i [path to reference index, if exists] \
   -o [output directory] \
   -s [which steps to skip: either or combination of characters {cudqm} ]

*Note: this should be a single-line command or the backslashes should be properly indented.*

The option ``-s`` specifies which steps of the pipeline you'd like to skip. Possible values are:

- c skip trimming adapters with cutadapt

- u skip UMI extraction

- d skip deduplication after alignment

- q skip quality check: FASTQC and MULTIQC

- m skip mapping

You may use any combination of these characters, e.g. use ``-s cudqm`` to skip all

This script will produce subfolders in the output directory, containing results of each step of the pipeline. The bam files will be generated in the **align_dedup** folder. 

In addition to performing the steps described above, it also evaluates the distribution of reads across the genome, according to gene classes {"rRNA" "mRNA" "tRNA" "snoRNA" "snRNA" "ncRNA"}. These statistics are kept in the **align_rna/rna_stats.txt** file. 

**!!NOTE!!** This example pipeline treats files as **singl-end** libraries. If you have paired-end reads, you should only supply the first read (\*_R1\* files) to fivepseq. 



Test data
---------

You may download the preprocessing test dataset for yeast `from here <http://data.pelechanolab.com/software/fivepseq/data/test_data_1.0b5/preprocess/>`_ to start with the following fastq to bam preprocessing pipeline. The archive contains subsetted fastq files from Saccharomyces Cerevisiae (first 100K reads), as well as genome fasta and annotation gff files. If you use another organism, you may refer to the following section to find needed genome files. 

The src folder in the test data contains a script named **preprocess_test.sh**, which contains the following command: 

.. code-block:: shell

 ./src/fivepseq_preprocess.sh \
     -f fastq \
     -o preprocess \
     -g genome/* \
     -a gff/* 
 
Navigate to the test data directory and give it a try by running ``./src/preprocess.sh``. It should generate the following folders in the preprocess directory: 

:fastqc_raw: FASTQC reports on the fastq files before preprocessing
:multiqc_raw: MULTIQC report on the fastq files before preprocessing
:fastq_trimmed: fastq files with standard illumina adapters :fastq_deumi: fastq files with UMI sequences extracted from the reads and placed in the read headers. Use ``-s ud`` if your data does not contain UMIs. 
:fastqc: FASTQC reports on the fastq files after trimming and UMI extraction
:multiqc: MULTIQC report on the fastq files after trimming and UMI extraction. 
:star_index: the reference genome index generated by STAR. Use ``-i star_index`` option in the future runs to use the already generated index. 
:star_align: alignment bam files 
:align_dedup: alignment bam files after UMI-based deduplication. These are the final files you should give as input to fivepseq. 
:filtered_gff: annotation files filtered by RNA content, used for overlapping with bam files to analyze RNA content (see below). 
:align_rna: bam files, filtered for the {"rRNA" "mRNA" "tRNA" "snoRNA" "snRNA" "ncRNA"} classes of RNA and the **rna.stats.txt** file that summarized the read counts for each RNA class. 
:bedgraph: two bedtraph files for 5'P endpoints mapping to fwd and rev strands for each bam file 
