Welcome to fivepseq readme!
==========================

Fivepseq is a software package for analysis of 5′ endpoints distribution in RNA sequencing datasets.
 

Installation
--------
Install dependencies:

To set up fivepseq, the following python packages need to be pre-installed manually using pip (if you don't have pip you may install it as described `here <ht\
tps://pip.pypa.io/en/stable/installing/>`_ ).

Paste the following lines into the shell terminal:

.. code-block:: shell

 pip install --upgrade numpy pysam cython
 pip install plastid


Clone the project from github:

.. code-block:: shell

 git clone https://github.com/lilit-nersisyan/fivepseq.git

Navigate into the fivepseq directory and install:

.. code-block:: shell

 python setup.py install

To check if fivepseq was installed correctly, type the following in the command line:

.. code-block:: shell

 fivepseq --version

This should display the currently installed version of fivepseq. To display commandline arguments you may type:

.. code-block:: shell

 fivepseq --help


In order to enable exporting vector and portable image files, you'll also need to have phantomjs installed as follows:

.. code-block:: shell

 conda install phantomjs selenium pillow




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

For UPPMAX users only
--------

- Install the latest stable version of fivepseq by
- - cd /proj/sllstore2017018/lilit/fivepseq_latest:
- - python setup.py install

*****************
Running fivepseq
*****************

Fivepseq requires the following files to run:

|    Aligned reads (**.bam**)
|    Alignment index (**.bai**)
|    Genomic sequence file (**.fasta** / **.fa**)
|    Genomic annotation file (**.gff/ .gtf**)

This section assumes that you already have these files. If not, please, refer to the section: **Preparing data**.

Fivepseq usage
----------------------------------

The ``fivepseq --help`` command will show fivepseq usage and will list all the arguments.

.. code-block:: shell

 usage: fivepseq -b B -g G -a A [optional arguments]

Required arguments
==================

.. code-block:: shell

 -b B   the full path one or many bam/sam files (many files should be provided with a pattern, **within double quotes**: e.g. ["your_bam_folder/*.bam"])
 -g G   the full path to the fa/fasta file
 -a A   the full path to the gtf/gff/gff3 file

**Note:**

- The indexed alignment files should be in the same directory as bam files, with the same name, with .bai extension added.

- Multiple bam files should be indicated with a pattern placed **within double quotes**: e.g. ["your_bam_folder/\*.bam"]


Commonly, you will run fivepseq by also providing the name of the **output folder** ('fivepseq' by default) and the **title** of your run (determined from bam path otherwise):

.. code-block:: shell

 fivepseq \
    -g <path_to_genome_fasta> \
    -a <path_to_annotation> \
    -b <path_to_bam_file(s) \
    -o <output_directory> \
    -t <title_of_the_run>

*Note: this is a single commandline, the backslashes are used to move to a new line for cozy representation: either copy-paste like this or use a single line without the backslashes.*

Optional arguments
==================

.. code-block:: shell

 --span SPAN

This specifies the number of bases to span around translation START and STOP positions. The default value is 100 bases.

.. code-block:: shell

 --conflicts {add,overwrite,alt_dir}

The conflict mode specifies how to deal with files/folders that already exist. You may choose either of the following options:

:add: (default). Only adds missing files to an existing directory.
:overwrite: Overwrites all the files in an existing output directory
:alt_dir: Creates an alternative directory by appending '+' suffix to an existing output directory

.. code-block:: shell

 --ignore-cache

When the annoatation gff/gtf file is read by fivepseq for the first time, it stores the transcript assembly object it generates in a pickle path, which is located in the same parent directory where your fivepseq output folder is, in the folder **fivepseq_cache**. If the annotation name stays the same, fivepseq will directly load this object in all further runs, instead of processing the annotation file. However, if the content of your annotation file has changed or you suspect that the previous pickle object might be truncated, you can use this option to ignore the cache and process the annotation file from scratch.

Additional arguments
====================

.. code-block:: shell

 -gf GF, -genefilter GF

If you are interested only in a specific set of transcripts, you may specify them with a text file containing newline-separated list of names you'd like to use. Note: only these genes will be used in all the calculations.

In a properly formatted gene filter file, the first line should specify the gene attribute in the gff file, and the rest of the lines should correspond to the actual values of that attribute in the gff file. For example, for the following two entries in the gff file:

::

 Isgdgene335649.+.ID=gene:YAL069W;biotype=protein_coding;gene_id=YAL069W;logic_name=sgd
 IsgdmRNA335649.+.ID=transcript:YAL069W_mRNA;Parent=gene:YAL069W;biotype=protein_coding;transcript_id=YAL069W_mRNA
 Isgdexon335649.+.Parent=transcript:YAL069W_mRNA;Name=YAL069W_mRNA-E1;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=YAL069W_mRNA-E1;rank=1
 IsgdCDS335649.+0ID=CDS:YAL069W_mRNA;Parent=transcript:YAL069W_mRNA;protein_id=YAL069W_mRNA

 Isgdgene538792.+.ID=gene:YAL068W-A;biotype=protein_coding;gene_id=YAL068W-A;logic_name=sgd
 IsgdmRNA538792.+.ID=transcript:YAL068W-A_mRNA;Parent=gene:YAL068W-A;biotype=protein_coding;transcript_id=YAL068W-A_mRNA
 Isgdexon538792.+.Parent=transcript:YAL068W-A_mRNA;Name=YAL068W-A_mRNA-E1;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=YAL068W-A_mRNA-E1;rank=1
 IsgdCDS538792.+0ID=CDS:YAL068W-A_mRNA;Parent=transcript:YAL068W-A_mRNA;protein_id=YAL068W-A_mRNA

The following gene filter file may be used, where the attribute is **gene_id**.

.. code-block:: shell

 gene_id
 YAL069W
 YAL068W-A

*In fact, the gene IDs you see in the file have the form* **gene:YAL069W**, *however, fivepseq tolerates if you just give the identifier following the colon.*

You may also specify the CDS ID, again omitting (or not) the text before colon:

.. code-block:: shell

 ID
 YAL069W_mRNA
 YAL068W-A_mRNA

You may also use other attributes in the CDS entry if you have alternative gff/gtf files. In case you specify the wrong attribute, the error message will tell you what attributes you may use. In case you use the wrong attribute values, the error message will list a few values of that attribute in the correct format.
fivepsq

.. code-block:: shell

 -gs GS, -genesets GS

If you'd like to compare fivepseq plots between different gene sets, you may provide this as an additional argument. In comparison to the gene filter option above, the gene set option works as follows. Fivepseq first performs regular calculations, including all the genes in the annotation file, and later generates additional reports to compare profiles based on the gene sets provided.

The file should be tab-delimited. With the first column indicating the genes, while the second one - the gene set names. The first line again indicates the attribute name and the geneset headding, while the rest of the lines contain gene - gene set mappings. The following example demonstrates one such file (the attribute choice is explained above):

.. code-block:: shell

   gene_id   GO:BP
   YLR116W   mRNA splicing, via spliceosome
   YDL070W   DNA repair
   YPR030W   fungal-type cell wall organization
   YJL158C   fungal-type cell wall organization
   YFL039C   DNA repair


Advanced arguments
==================

Fivepseq reduces noise by detecting 5' counts that our outliers in the background count distribution. The latter is well approximated with Poisson distribution, computed based on the count distribution mean. Counts for which the probability of falling into this distribution is less than a certain threshold (0 by default) are considered as outliers. These outliers are usually down-scaled or down-sampled to the most extreme distribution count possible. However, you may modify this by either of the two options below:

.. code-block:: shell

    --ds DS, --downsample DS

With this option you can omit the distribution-modeled down-sampling described above and specify a constant value instead. Counts exceeding this threshold will be down-scaled to it.

.. code-block:: shell

    --op OP, --outlier_probability OP

With this option you may change the default probability threshold of Poisson distribution. You may increase it from 0, to be more harsh in allowing high count values.
You can also use this option to turn off down-sampling altogether, by setting the probability threshold to -1.


Obsolete
==================

**--loci-file**

*default: None*

This option requires a file with coordinates of the loci (e.g. RBP binding coordinates), relative to which, the user wants to generate line-charts.
The file should be tab-separated, with the following structure:
Columns: chr->str->start->end->symbol
Rows: the chromosome name, strand (+ or -), start and end coordinates and the name of the RBP (or the locus).

This feature is in beta-, so line charts of reads relative to all the loci combined will be plotted underneath the main html file.
Four different plots will correspond to reads located in (1) 3UTR and CDS regions, (2) only 3UTR, (3) only 5UTR and (4) only CDS.


Have fun!
