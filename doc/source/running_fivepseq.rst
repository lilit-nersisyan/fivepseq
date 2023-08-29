.. _running_fivepseq:

*****************
Running fivepseq
*****************

Fivepseq requires the following files to run: 

|    Aligned reads (**.bam**)
|    Alignment index (**.bai**)
|    Genomic sequence file (**.fasta** / **.fa**)
|    Genomic annotation file (**.gff/ .gtf**)

This section assumes that you already have these files. If not, please, refer to the section: **Preparing the data**. 

Fivepseq usage
----------------------------------

The ``fivepseq --help`` command will show fivepseq usage and will list all the arguments. 

.. code-block:: shell

 usage: fivepseq -b B -g G -a A [optional arguments]

Required arguments
==================

.. code-block:: shell

 -b B   full path to one or more bam/sam files (many files should be provided with a pattern, **within double quotes**: e.g. ["your_bam_folder/*.bam"])
 -g G   full path to the fa/fasta file
 -a A   full path to the gtf/gff/gff3 file

**Note:** 

- The indexed alignment files (.bam.bai) should be in the same directory as the bam files, with the same name and .bai extension. 

- Multiple bam files should be indicated with a pattern placed **within double quotes**: e.g. ["your_bam_folder/\*.bam"]


Commonly, you will run fivepseq by also providing the name of the **output folder** ('fivepseq' by default) and the **title** of your run (if not provided, it will be determined from bam path): 

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

The conflicts mode specifies how to deal with files/folders that already exist. You may choose either of the following options:

:add: (default). Only adds missing files to an existing directory and skips generating the ones that already exist. 
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
 
.. code-block:: shell

 -gs GS, -genesets GS

If you'd like to compare fivepseq plots between different gene sets, you may provide this as an additional argument. In comparison to the gene filter option above, the gene set option works as follows. Fivepseq first performs regular calculations, including all the genes in the annotation file, and later generates additional reports to compare profiles based on the gene sets provided.  

The file should be tab-delimited. With the first column indicating the genes, while the second one - the gene set names. The first line again indicates the attribute name and the geneset heading, while the rest of the lines contain gene - gene set mappings. The following example demonstrates one such file (the attribute choice is explained above): 

.. code-block:: shell

   gene_id   GO:BP
   YLR116W   mRNA splicing, via spliceosome
   YDL070W   DNA repair
   YPR030W   fungal-type cell wall organization
   YJL158C   fungal-type cell wall organization
   YFL039C   DNA repair


*Note: the first line is required.* 

Advanced arguments 
==================
**Options for codon-specific counts**

The gene boundaries are usually masked by 20 nt (default) when counting the distribution of 5'P counts relative to codons and codon-motifs. You may change the mask size by setting it in the [3:50) range.

.. code-block:: shell

    -ms [mask size],  -codon-mask-size [mask size]

Or you may turn masking off altogether, by passing in the following option:

.. code-block:: shell

    --no-mask

After analyzing di/tri-codon/peptide specific counts *fivepseq* chooses top 50 motifs that have highest pausing at the position which corresponds to ribosome protected fragments when the last codon of the motif is in the A site, which is -14nt for di-peptides(codons) and -11nt for tri-peptides(codons). You may change these positions with the following options:

.. code-block:: shell

    -dipeptide-pos

Pass a number in the range [-27:6) to this option to tell fivepseq counts in which position {-27:6} from A site should be ordered to output top stalled dipeptides.

.. code-block:: shell

    -tripeptide-pos

Pass a number in the range [-24:9) to this option to tell fivepseq counts in which position {-24:9} from A site should be ordered to output top stalled tripeptides.

**Noise removal**

Fivepseq reduces noise by detecting 5' counts that our outliers in the background count distribution. The latter is well approximated with Poisson distribution, computed based on the count distribution mean. Counts for which the probability of falling into this distribution is less than a certain threshold (0 by default) are considered as outliers. These outliers are usually down-scaled or down-sampled to the most extreme distribution count possible. However, you may modify this by either of the two options below:

.. code-block:: shell

    --ds DS, --downsample DS

With this option you can omit the distribution-modeled down-sampling described above and specify a constant value instead. Counts exceeding this threshold will be down-scaled to it.  

.. code-block:: shell

    --op OP, --outlier_probability OP

With this option you may change the default probability threshold of Poisson distribution. You may increase it from 0, to be more harsh in allowing high count values. 
You can also use this option to turn off down-sampling altogether, by setting the probability threshold to -1. 

**Transcript features**

By default fivepseq filters genes from the annotation file that have mRNA or protein_coding tags. You may change it to your desired type available under 'type' or 'biotype' attributes in your annotation.

.. code-block:: shell

    -transcript-type


**Subset**

Fivepseq will only perform analysis for the first genes, the number being equal to the argument passed to this option.

.. code-block:: shell

    -subset

**Plotting options**


.. code-block:: shell

    -triangle_threshold

Count threshold (0:10000) for points in the triangle plot (points with lower counts will not be plotted)


Running fivepseq on test data
-----------------------------
To run fivepseq on test alignment files, download the dataset `from here <http://data.pelechanolab.com/software/fivepseq/data/test_data_1.0b5/yeat_eIF5A/>`_. 

The test data folder also contains source files for running fivepseq. There is a script for a quick test (**src/fivepseq_quick.sh** in the test_data folder), where only a small subset of transcripts is chosen. This test does not produce reasonable plots, but if it runs successfully you will know that everything is set up properly. Navigate to the test dataset directory (yeast_eIF5A) and run the **./src/fivepseq_quick.sh** script, or copy paste this code to the command line:
    
.. code-block:: shell

 fivepseq \
    -g genome/* \
    -a gff/* \
    -b "preprocess/align_dedup/*.bam" \
    -o fivepseq_quick \
    -t yeast_quick_test \
    -gf genesets/test_set.txt

This should produce a folder named **fivepseq_quick**, with the following file structure: 

.. code-block:: shell

    fivepseq_quick
    ├── fivepseq_counts
    │   ├── 1_eIF5A_mut_R1
    │       ├── protein_coding
    │           ├── text files explained in the next chapters
    │   ├── 2_WT_R3
    │       ├── protein_coding
    │           ├── text files explained in the next chapters
    │   ├── 3_eIF5A_mut_R1
    │       ├── protein_coding
    │           ├── text files explained in the next chapters
    ├── fivepseq_plots
    │   ├── main
    │       ├── png
    │       ├── svg
    │       ├── yeast_quick_test_combined.html
    │       ├── yeast_quick_test_main.html
    │   ├── supplement
    │       ├── png
    │       ├── svg
    │       ├── yeast_quick_test_amino_acid_linecharts.html
    │       ├── yeast_quick_test_codon_heatmaps.html
    │       ├── yeast_quick_test_codon_linecharts.html
    │       ├── yeast_quick_test_dicodon_linecharts.html
    │       ├── yeast_quick_test_dipeptide_linecharts.html
    │       ├── yeast_quick_test_tricodon_linecharts.html
    │       ├── yeast_quick_test_tripeptide_linecharts.html
    │   ├── comparison
    │       ├── yeast_quick_test_differential_heatmaps.html
    │   ├── yeast_quick_test.html   # <<-- this is the main entry point for all the reports 
    ├── log
    │   ├── fivepseq.log
    ├── count_summary.txt

The main result of the run is the file **fivepseq_quick/fivepseq_plots/yeast_quick_test.html**. This file has hyperlinks to all the other report files. The files generated here will be explained in detail in the **Interpreting the output** section. 

