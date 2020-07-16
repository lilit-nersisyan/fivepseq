Welcome to fivepseq readme!
==========================

Fivepseq is a software package for analysis of 5′ endpoints distribution in RNA sequencing datasets.
 

Homepage
------------
The homage is hosted at Pelechano lab website at http://pelechanolab.com/software/fivepseq/.

User guide
------------
Below is a quick manual to get you started.
For detailed instructions and explanations on fivepseq output, please visit the user guide at: https://fivepseq.readthedocs.io/en/latest/.

Installation
------------
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


Running fivepseq
==================

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

Additional arguments
==================

Type fivepseq --help to see the list of additional arguments. For a detailed description of available arguments, see the User guide at: https://fivepseq.readthedocs.io/en/latest/.



Preprocessing from FASTQ files
----------------------------------------
Fastq files need to be preprocessed and aligned to the reference genome before proceeding to fivepseq downstream analysis. Preprocessing proceeds with the following steps:

- quality checks (with FASTQC and MULTIQC),
- adapter and quality based trimming,
- UMI extraction (if the library was generated with UMIs),
- mapping to reference
- read deduplication (if the library was generated with UMIs),
- bedgraph generation to view 5'P count distribution in genome viewers

An example of pre\-processing pipeline can be found in the preprocess_scripts directory

In order to run this pipeline, you need to have access to common bioinformatics software such as `STAR <https://github.com/alexdobin/STAR>`_, `UMI-tools <https://github.com/CGATOxford/UMI-tools>`_, `bedtools <https://bedtools.readthedocs.io/en/latest/>`_, `Samtools <http://www.htslib.org/>`_, `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_, `MultiQC <https://multiqc.info/>`_ and `cutadapt <https://github.com/marcelm/cutadapt>`_.

To use it, navigate to the directory where the script is located and use the following command in the prompt:

.. code-block:: shell

 ./fivepseq_preprocess.sh -f [path to directory containing fastq files] -g [path to genome fasta] -a [path to annotation gff/gtf] -i [path to reference index, if exists] -o [output directory] -s [which steps to skip: either or combination of characters {cudqm} ]

The option ``-s`` specifies which steps of the pipeline you'd like to skip. Possible values are:

- c skip trimming adapters with cutadapt

- u skip UMI extraction

- d skip deduplication after alignment

- q skip quality initial check: FASTQC and MULTIQC

- p skip post-processing quality check: FASTQC and MULTIQC

- m skip mapping

- d skip deduplication

You may use any combination of these characters, e.g. use ``-s cudqm`` to skip all

This script will produce sub-folders in the output directory, containing results of each step of the pipeline. The bam files will be generated in the **align_dedup** folder.

In the  In addition to performing the steps described above, it also evaluates the distribution of reads across the genome, according to gene classes {"rRNA" "mRNA" "tRNA" "snoRNA" "snRNA" "ncRNA"}. These statistics are kept in the **align_rna/rna_stats.txt** file.

**!!NOTE!!** This example pipeline treats files as **singl-end** libraries. If you have paired-end reads, you should only supply the first read (\*_R1\* files) to fivepseq.


For UPPMAX users only
------------------------

- Install the latest stable version of fivepseq by
- - cd /proj/sllstore2017018/lilit/fivepseq_latest:
- - python setup.py install

Have fun!
