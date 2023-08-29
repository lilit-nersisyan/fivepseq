.. fivepseq documentation master file, created by
   sphinx-quickstart on Tue Dec  3 18:56:05 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Fivepseq User Guide
====================

Fivepseq is a software package for analysis of 5' endpoint distribution in RNA sequencing datasets. This is particularly useful for techniques that capture 5' monophosphorylated RNAs, such as 5PSeq, PARE-seq or GMUC. It may also be useful for ribosome profiling datasets and alike. 

The main workflow of fivepseq is intended for downstream analysis of alignment files to describe the distribution of 5' endpoints of reads relative to translation start and stop sites, as well as relative to amino acids or codons. It also computes frame preference of 5' endpoint distribution and captures periodicity patterns. 

The output of fivepseq is in the form of html report files, as well as text files and image files in portable network graphic and vector graphics formats. 

Homepage
_________

For more information, demo cases and citation, visit `Fivepseq homepage <http://pelechanolab.com/software/fivepseq/>`_

==================

.. toctree::
   :maxdepth: 2

   installation.rst
   preparing_data.rst
   running_fivepseq.rst
   interpreting_output.rst
   count_files.rst
