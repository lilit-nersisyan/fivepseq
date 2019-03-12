Welcome to fivepseq readme!  
--------
This program reports 5'-footprint properties from 5Pseq reads in a given alignment file. 
 
There are to modes to run fivepseq: 
count - generates count vectors and meta-vectors from the given alignment file and stores them in the output directory as text files
plot - generates plots based on input count data folder and stores those to an html file

Usage
--------
fivepseq call should start with fivepseq path followed by required arguments common to both count and plot commands, to be followed by either count or plot specification with their specific arguments. 

python path_to_your_fivepseq (/proj/sllstore2017018/lilit/fivepseq_v0.1.1/fivepseq) \
-g path_to_your_genome_fasta \
-a path_to_your_annotation_gff (sorry no gtf for now) \
[--conflicts conflict_mode]

optional argument: --conflicts

The conflict mode specifies how to deal with files/folders that alraedy exist. There are three options you may choose from: 

- add (default) - only adds missing files to existing directory 

- overwrite - overwrites all the files in the existing (in case) output directory 

- alt_dir - uses alternative directory by appending '+' suffix to existing (in case) output directory



Usage - count
-----

To call count, simply specify the count command after the common syntax described above:

python path_to_your_fivepseq (/proj/sllstore2017018/lilit/fivepseq_v0.1.1/fivepseq) \
-g path_to_your_genome_fasta \
-a path_to_your_annotation_gff (sorry no gtf for now) \
count \
-o path_to_output_directory


Usage - plot
-----
To call count, simply specify the count command after the common syntax described above:

python path_to_your_fivepseq (/proj/sllstore2017018/lilit/fivepseq_v0.1.1/fivepseq) \
-g path_to_your_genome_fasta \
-a path_to_your_annotation_gff (sorry no gtf for now) \
plot \
-sd or -md path_to_count_folder(s)
-o path_to_output_directory
-t title_of_html_file

Note that plot function can take as input a single (-sd) or multiple (-md) count directories. 

For example, if you have single count folder then you should specify: 
-sd my_one_and_only_count_folder
-t my_one_and_only_sample

If you have multiple count folders, in a parent_directory and you'd like to generate plots for several of them in one html file, you should give those with -md option. You can also specify the cound directories with a pattern: 

-md "my_parent_count_directory/all_folders_starting_with_S_cer*"
-t S_cer

!!! Note the brackets "" following the -md command: don't miss them when specifying multiple directories 


Have fun! 
