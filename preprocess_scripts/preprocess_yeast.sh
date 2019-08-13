#!/bin/sh

#SBATCH -A snic2018-8-292
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 17:00:00
#SBATCH -J preprocesing_fivepseq_human
#SBATCH -o log/preprocesing_fivepseq_human.log

echo "################################" 
echo "###         INPUT            ###" 
echo "################################" 
echo ""

while getopts f:o:s:g:a:i:d option
do
case "${option}"
in
f) FASTQ_DIR=${OPTARG};;
o) OUT_DIR=${OPTARG};;
s) SKIP=${OPTARG};;
g) GENOME=$OPTARG;;
a) ANNOT=$OPTARG;;
i) INDEX=$OPTARG;;
d) 
esac
done


fastq_dir=$FASTQ_DIR
out_dir=$OUT_DIR
skip=$SKIP

if [ -z $GENOME ]; then
    echo "ERROR: no genome fasta file supplied"
    exit 1
fi 

if [ ! -f $GENOME ]; then 
    echo "ERROR: the genome fasta file $GENOME does not exist"
    exit 1
fi 

if [ -z $ANNOT ]; then
    echo "ERROR: no annotation gff file supplied"
    exit 1
fi 

if [ ! -f $ANNOT ]; then 
    echo "ERROR: the genome fasta file $ANNOT does not exist"
    exit 1
fi 


echo "Supplied arguments:"
echo "FASTQ DIRECTORY: $fastq_dir"
echo "OUTPUT DIRECTORY: $out_dir"
echo "SKIP: $skip"
echo "GENOME FASTA: $GENOME"
echo "ANNOTATION: $ANNOT"
echo "INDEX: $INDEX"

if grep -q 'c' <<<"$skip"; then
    echo ""
    echo "Adapter trimming will be skipped"
fi

if grep -q 'u' <<<"$skip"; then
    echo ""
    echo "UMI extraction will be skipped"
fi

if grep -q 'd' <<<"$skip"; then
    echo ""
    echo "Read deduplication will be skipped"
fi

if grep -q 'q' <<<"$skip"; then
    echo ""
    echo "QC will be skipped"
fi

if grep -q 'm' <<<"$skip"; then
    echo ""
    echo "Mapping will be skipped"
fi

if [ ! -d $out_dir ]; then 
    mkdir $out_dir
else
    echo ""
    echo "WARNING: the contents of the directory $out_dir may be overwritten." 
    echo ""
fi


if [ -d $fastq_dir ]
then     
    n_files=`ls $fastq_dir/*.gz | wc -l`
    if [ $n_files -eq 0 ]; then
	echo "ERROR: No compressed fastq files found. Exiting."
	exit 1
    else
    
	echo "The following fastq gzipped files will be processed:"
	for f in $fastq_dir/*.gz
	do
	    echo "$f"
	done
   fi
else 
    echo "ERROR: the provided fastq file $fastq_in does not exist" 
    exit 1
fi 

echo ""


if grep -q 'q' <<<"$skip"; then
    echo ""
    echo "################################" 
    echo "###        SKIP: QC          ###" 
    echo "################################" 
    echo ""

else

    echo ""
    echo "################################" 
    echo "###         FASTQC           ###" 
    echo "################################" 
    echo ""

    fastqc_dir=$out_dir/fastqc_raw


    echo "INFO: fastqc input" 
    echo "  fastq: $fastq_dir"
    echo "  output: $fastqc_dir"
    echo ""

    if [ ! -d $fastqc_dir ]; then 
	mkdir $fastqc_dir
    else
	echo ""
	echo "WARNING: the contents of the directory $fastqc_dir may be overwritten." 
	echo ""
    fi

    module load bioinfo-tools FastQC
    wait

    fastqcdir=`which fastqc`
    if [ -z $fastqcdir ]; then
	echo "Could not load FastQC"
    else
	echo ""
    fi


    for f in $fastq_dir/*.gz
    do
	echo "$f: fastqc"
	fastqc -o $fastqc_dir $f
	wait
    done
    echo "FastQC finished for all fastq files in $fastq_dir. Output written to $out_dir"


    echo ""
    echo "################################" 
    echo "###         MULTIQC           ###" 
    echo "################################" 
    echo ""

    multiqc_dir=$out_dir/multiqc_raw

    echo "INFO: fastqc input" 
    echo "  fastqc: $fastqc_dir"
    echo "  output: $multiqc_dir"
    echo ""


    if [ ! -d $multiqc_dir ]; then 
	mkdir $multiqc_dir
    else
	echo ""
	echo "WARNING: the contents of the directory $multiqc_dir may be overwritten." 
	echo ""
    fi

    module load bioinfo-tools MultiQC
    wait
    multiqcdir=`which multiqc`
    if [ -z $multiqcdir ]; then
	echo "Could not load MULTIQC"
    else
	echo ""
    fi

    multiqc -d -o $multiqc_dir $fastqc_dir


    echo "MULTIQC finished for all fastqc files in $fastqc_dir"

    module unload MultiQC
fi
if grep -q 'c' <<<"$skip"; then
    echo ""
    echo "################################" 
    echo "###        SKIP: CUTADAPT    ###" 
    echo "################################" 
    echo ""

else

    echo "################################" 
    echo "###         CUTADAPT         ###" 
    echo "################################" 
    echo ""

    module load bioinfo-tools
    wait
    module load cutadapt
    wait

    moduledir=`which cutadapt`
    if [ -z $moduledir ]; then
	echo "Could not load cutadapt"
    else
	echo ""
    fi



    cutadapt_dir=$out_dir'/fastq_trimmed'
    echo "INFO: cutadapt input" 
    echo "  fastq: $fastq_dir"
    echo "  output: $cutadapt_dir"
    echo "  options: --minimum-length 28; -e 0.2; -o 9; --nextseq-trim 20"
    echo ""

    if [ ! -d $cutadapt_dir ]; then 
	mkdir $cutadapt_dir
    else
	echo "WARNING: the contents of the directory $cutadapt_dir may be overwritten." 
    fi

    echo ""

    for f in $fastq_dir/*.gz
    do
	fname=${f##*/}
	echo ""
	echo "$f: trimming"
	cutadapt -a AGATCGGAAGAGCAC --minimum-length 28 -e 0.2 -o 9 --nextseq-trim=20 -o $cutadapt_dir/$fname $f
	wait
	success=$?

	if [ $success -eq 0 ]; then 
	    echo "$f: successfully trimmed"
	    echo ""
	else
	    echo "Cutadapt finished with non-zero exit status $success for file $f"
	    exit $success
	fi
    done

    echo "Successfully finished adapter trimming for all files"
    echo ""
    fastq_dir=$cutadapt_dir
fi


if grep -q 'u' <<<"$skip"; then
    echo ""
    echo ""
    echo "################################" 
    echo "###     SKIP: UMI EXTRACT    ###" 
    echo "################################" 
    echo ""

else


    echo "################################" 
    echo "###         UMI EXTRACT      ###" 
    echo "################################" 
    echo ""
    

    module load bioinfo-tools
    wait
    module load umi_tools
    wait

    umidir=`which umi_tools`
    if [ -z $umidir ]; then
	echo "Could not load umi_tools"
    else
	echo ""
    fi

    echo ""

    deumi_dir=$out_dir'/fastq_deumi'

    echo "INFO: umitools input" 
    echo "  fastq: $fastq_dir"
    echo "  output: $deumi_dir"
    echo "  options: --bc-pattern NNNNNNNN"
    echo ""


    if [ ! -d $deumi_dir ]; then 
	mkdir $deumi_dir
    else
	echo "WARNING: the contents of the directory $deumi_dir may be overwritten." 
    fi

    echo ""
    

    for f in $fastq_dir/*
    do
	fname=${f##*/}
	fname=${fname%.*}
	echo ""
	echo "$f: extracting UMI"
	gunzip -c $f | umi_tools extract --bc-pattern NNNNNNNN --log $deumi_dir/$fname.umi.extract.log | gzip -c > $deumi_dir/$fname.gz
	wait
	success=$?

	if [ $success -eq 0 ]; then 
	    if [ ! -z $deumi_dir/$fname.gz ]; then 
		echo "$f: successfully extracted UMI"
		echo ""
	    else
		echo "UMI extraction finished with non-zero exit status $success for file $f"
		exit 2
	    fi		
	else
	    echo "UMI extraction finished with non-zero exit status $success for file $f"
	    exit $success
	fi
    done
    echo "UMI extraction finished for all fastq files in $cutadapt_dir. Output written to $deumi_dir"
    fastq_dir=$deumi_dir
fi 

if grep -q 'q' <<<"$skip"; then
    echo ""
    echo "################################" 
    echo "###      SKIP: (post) QC     ###" 
    echo "################################" 
    echo ""

else

    echo ""
    echo "################################" 
    echo "###     (post)  FASTQC       ###" 
    echo "################################" 
    echo ""

    fastqc_dir=$out_dir/fastqc


    echo "INFO: fastqc input" 
    echo "  fastq: $fastq_dir"
    echo "  output: $fastqc_dir"
    echo ""

    if [ ! -d $fastqc_dir ]; then 
	mkdir $fastqc_dir
    else
	echo ""
	echo "WARNING: the contents of the directory $fastqc_dir may be overwritten." 
	echo ""
    fi

    module load bioinfo-tools FastQC
    wait

    fastqcdir=`which fastqc`
    if [ -z $fastqcdir ]; then
	echo "Could not load FastQC"
    else
	echo ""
    fi


    for f in $fastq_dir/*.gz
    do
	echo "$f: fastqc"
	fastqc -o $fastqc_dir $f
	wait
    done
    echo "FastQC finished for all fastq files in $fastq_dir. Output written to $out_dir"

    echo ""
    echo "################################" 
    echo "###      (post) MULTIQC      ###" 
    echo "################################" 
    echo ""

    multiqc_dir=$out_dir/multiqc

    echo "INFO: fastqc input" 
    echo "  fastqc: $fastqc_dir"
    echo "  output: $multiqc_dir"
    echo ""


    if [ ! -d $multiqc_dir ]; then 
	mkdir $multiqc_dir
    else
	echo ""
	echo "WARNING: the contents of the directory $multiqc_dir may be overwritten." 
	echo ""
    fi

    module load bioinfo-tools MultiQC
    wait
    multiqcdir=`which multiqc`
    if [ -z $multiqcdir ]; then
	echo "Could not load MULTIQC"
    else
	echo ""
    fi

    multiqc -d -o $multiqc_dir $fastqc_dir


    echo "MULTIQC finished for all fastqc files in $fastqc_dir"

    module unload MultiQC
fi
if [ -f $INDEX/SAindex ]; then
    echo ""
    echo "################################" 
    echo "### SKIP: INDEX GENERATION   ###" 
    echo "################################" 
    echo ""

    echo "Supplied index directory $INDEX will be used as reference"
    genome_ind=$INDEX    
    echo ""
else
    echo "INFO: No proper STAR index was supplied. Will proceed with index generation."
    echo ""
    echo "######################################" 
    echo "###     GENERATE GENOME INDEX      ###" 
    echo "######################################" 
    echo ""

    if [ -f "$GENOME" ]; then 
	echo "Genome fasta file: $GENOME"
    else
	echo "ERROR: Genome fasta file was not supplied or $GENOME does not exist"    
	exit 1 
    fi
    
    if [ -f "$ANNOT" ]; then 
	echo "Genome annotation file: $ANNOT"
    else
	echo "ERROR: Genome annotation file was not supplied or $ANNOT does not exist"    
	exit 1 
    fi

    module load bioinfo-tools star
    wait

    echo "Loading STAR"
    stardir=`which star`
    if [ -z $stardir ]; then
	echo "ERROR: Could not load star"
	exit 3
    else
	echo ""
    fi

    genome_ind=$out_dir/star_index
    genome_fasta=$GENOME
    genome_gff=$ANNOT


    echo "INFO: star generate genome input" 
    echo "  genome fasta: $genome_fasta"
    echo "  genome gff: $genome_gff"
    echo "  output: $genome_ind"
    echo "  options: --sjdbGTFtagExonParentTranscript Parent (GFF annotation)"
    echo ""

    
    if [ -d $genome_ind ]; then
	echo "Index directory $genome_ind already exists. Deleting content."
	rm $genome_ind/*
	wait
    else
	echo "Creating directory $genome_ind"
	mkdir $genome_ind
    fi

    if [ ! -d $genome_ind  ]; then
	echo "ERROR: Could not create directory $genome_ind. Exiting from star.generate.genome.sh"
	exit 2
    fi

    echo "STAR: generating genome index"

    star \
	--runMode genomeGenerate \
	--genomeDir $genome_ind \
	--genomeFastaFiles $genome_fasta \
	--sjdbGTFfile $genome_gff \
	--sjdbGTFtagExonParentTranscript Parent 

    wait
    if [ -f $genome_ind'/SAindex' ]; then
	echo "Successfully generated STAR genome in directory: $genome_ind"
    
    else
	echo "ERROR: Problem generating STAR index: check the log files."
	exit 4
    fi
fi

if grep -q 'm' <<<"$skip"; then
    echo ""
    echo "################################" 
    echo "###      SKIP: MAPPING       ###" 
    echo "################################" 
    echo ""

    align_dir=$out_dir/star_align

    if [ ! -d $align_dir ]; then
	echo "ERROR: the skip option includes mapping, but the directory $align_dir does not exist"
	exit 6
    fi

else

    echo "######################################" 
    echo "###             MAPPING            ###" 
    echo "######################################" 

    module load bioinfo-tools star
    stardir=`which star`
    if [ -z $stardir ]; then
	echo "ERROR: Could not load star"
	exit 3
    else
	echo ""
    fi

    module load samtools
    wait

    samtoolsdir=`which samtools`
    if [ -z $samtoolsdir ]; then
	echo "Could not load samtools"
	exit 3
    else
	echo ""
    fi

    fastq_dir=$fastq_dir
    genome_ind=$genome_ind
    align_dir=$out_dir/star_align


    echo "INFO: star alignment input" 
    echo "  fastq: $fastq_dir"
    echo "  genome index: $genome_ind"
    echo "  output: $align_dir"
    echo "  options: --alignEndsType Extend5pOfRead1; --outFilterMatchNminOverLread 0.9; --outFilterMultimapNmax 2; --limitBAMsortRAM 100000000000; alignIntronMax 2500"
    echo "  comment: primary alignments are chosen separately; unmapped reads are kept in a separate file"
    echo ""

    if [ ! -d $align_dir ];then
	echo "Creating output directory $align_dir"
	mkdir $align_dir
	wait
    else
	echo "WARNING:: The output directory $align_dir already exists. It's content will be overwritten.\n"
    fi


    for f in $fastq_dir/*.gz
    do
	fname=${f##*/}
	fname=${fname%.*}
	echo ""
	echo "$f: STAR alignment"
	out_prefix=$align_dir/$fname
	star \
            --readFilesIn $f --readFilesCommand gunzip -c \
            --genomeDir $genome_ind \
            --outFileNamePrefix $out_prefix \
            --outSAMtype BAM SortedByCoordinate \
            --alignEndsType Extend5pOfRead1 \
            --outFilterMatchNminOverLread 0.9 \
            --outFilterMultimapNmax 3 \
	    --alignIntronMax 2500 \
            --outReadsUnmapped Fastx \
            --limitBAMsortRAM 1000000000


	wait

	success=$?

	if [ $success -eq 0 ]; then 
	    echo ""
	else
	    echo "STAR alignment finished with non-zero exit status $success for file $fname"
	    exit $success
	fi

	echo "$f: STAR alignment done. Choosing primary alignments"
	samtools view -F 0x100 -b $out_prefix'Aligned.sortedByCoord.out.bam' > $out_prefix'Aligned.primary.out.bam'
	wait

	# samtools indexing

	if [ ! -f $out_prefix'Aligned.primary.out.bam' ]; then
	    echo "ERROR: Problem occured during STAR alignment for file $fname"
	    exit 4
	else
	    echo "$out_prefix'Aligned.primary.out.bam': indexing"
	    samtools index $out_prefix'Aligned.primary.out.bam'
	    wait
	    
	    success=$?
	    if [ $success -eq 0 ]; then 
		echo "$out_prefix'Aligned.primary.out.bam': successfully indexed"
		
	    else
		echo "Samtools indexing finished with non-zero exit status $success for file $fname"
		exit $success
	    fi

	fi
	
    done

    echo "INFO: alignment finished for all files in the directory $fastq_dir"
fi

if grep -q 'd' <<<"$skip"; then
    echo ""
    echo "################################" 
    echo "###        SKIP: DEDUP    ###" 
    echo "################################" 
    echo ""
    dedup_dir=$out_dir/align_dedup
    if [ ! -d $dedup_dir ]; then 
	echo "WARNING: deduplication directory not found. Skipping deduplication, using alignment files instead"
	dedup_dir=$align_dir
    fi
else

    echo "################################" 
    echo "###      DEDUPLICATION       ###" 
    echo "################################" 

    module load bioinfo-tools umi_tools
    wait


    umidir=`which umi_tools`
    if [ -z $umidir ]; then
	echo "Could not load umi_tools"
    else
	echo ""
    fi

    module load samtools
    wait

    samtoolsdir=`which samtools`
    if [ -z $samtoolsdir ]; then
	echo "Could not load samtools"
    else
	echo ""
    fi


    bam_dir=$align_dir
    dedup_dir=$out_dir/align_dedup

    echo "INFO: umitools deduplicate input" 
    echo "  bam: $align_dir"
    echo "  output: $dedup_dir"
    echo ""

    if [ ! -d $dedup_dir ];then
	echo "Creating output directory $dedup_dir"
	mkdir $dedup_dir
	wait
    else
	echo "WARNING:: The output directory $dedup_dir already exists. It's content will be overwritten.\n"
    fi


    for f in $bam_dir/*primary.out.bam
    do
	fname=${f##*/}
	fname=${fname%.*}
	echo ""
	echo "$fname: UMI deduplication"
	umi_tools dedup -I $f -S $dedup_dir/$fname"_dedup.bam" --log $dedup_dir/$fname".umi.dedup.log"
	wait

	echo "$dedup_dir/$fname'_dedup.bam': indexing"
	samtools index $dedup_dir/$fname"_dedup.bam"
	wait
	echo "$dedup_dir/$fname'_dedup.bam': successfully indexed"
    done
    echo "Deduplication finished for all bam files in $bam_dir. Output written to $dedup_dir"
fi


echo "################################" 
echo "###       RNA CONTENT        ###" 
echo "################################" 

module load bioinfo-tools BEDTools
wait


bedtoolsdir=`which bedtools`
if [ -z $bedtoolsdir ]; then
    echo "Could not load bedtools"
else
    echo ""
fi

module load bioinfo-tools samtools
wait
samtoolsdir=`which samtools`
if [ -z $samtoolsdir ]; then
    echo "Could not load samtools"
else
    echo ""
fi


gff_dir=$out_dir/filtered_gff
if [ ! -d $gff_dir ];then
    echo "Creating target gff directory $gff_dir"
    mkdir $gff_dir
    wait
else
    echo "WARNING:: The output directory $gff_dir already exists. It's content will be overwritten.\n"
fi


align_rna=$out_dir'/align_rna'
if [ ! -d $align_rna ];then
    echo "Creating classified alignment directory $align_rna"
    mkdir $align_rna
    wait
else
    echo "WARNING:: The output directory $align_rna already exists. It's content will be overwritten.\n"
fi

stats_file=$align_rna/rna.stats.txt
if [ -f $stats_file ]; then
    rm $stats_file
fi

echo "INFO: genome position stats input" 
echo "      gff: $ANNOT"
echo "      bam: $dedup_dir"
echo "  gff_dir: $gff_dir"
echo "   output: $align_rna"
echo "    stats: $stats_file"
echo ""


echo "Subsetting source gff $ANNOT"

rna_types=( "rRNA" "mRNA" "tRNA" "snoRNA" "snRNA" "ncRNA" )
for rna_type in ${rna_types[@]}
do 
    echo $rna_type
    cat $ANNOT | awk -v rna_type=$rna_type '$3 == rna_type { print $0 }' > $gff_dir'/'$rna_type'.gff3'

done 

echo "Intersecting alignment files"
for bam in $dedup_dir/*primary*.bam
do
    bamname=${bam##*/}
    bamname=${bamname%.*}
    echo $bamname
    
    for rna_type in ${rna_types[@]}
    do
	bam_rna=$align_rna/$bamname'_'$rna_type'.bam'
	bedtools intersect -s -abam $bam -b $gff_dir'/'$rna_type'.gff3' -wa > $bam_rna
	wait
	num_reads=$(samtools view $bam_rna | wc -l)
	echo "$bamname|$rna_type|$num_reads" >> $stats_file
    done
done



echo "       SUCCESS           "
