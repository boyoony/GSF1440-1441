#!/bin/bash

# GSF1440-1441

# 205 L1 N2 and 206 L1 adr-2(-) 
#GSF1440-WW205-2_S5_R1_001_trimmed.fq.gz
#GSF1440-WW205-3_S9_R1_001_trimmed.fq.gz
#GSF1440-WW206-2_S6_R1_001_trimmed.fq.gz
#GSF1440-WW206-3_S10_R1_001_trimmed.fq.gz

#GSF1441-WW205-1_S5_R1_001_trimmed.fq.gz
#GSF1441-WW206-1_S6_R1_001_trimmed.fq.gz


# 20210927_GSF1440_1441
Practicing RNA-seq analysis using L1 RNA-seq samples from our lab 


# Activate RNA-seq environment

conda activate rnaseq 


# Unload Perl

module unload perl


# Variables

ASSEMBLY='ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic.fa.gz'
ANNOTATION='ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.canonical_geneset.gtf.gz'


# Go to your slate account

cd /N/slate/by6


# Make a new folder 

mkdir GSF1440_1441


# Change directory to project space: 

cd /N/project/HundleyLab


# Secure copy and paste files from HundleyLab project space 

scp GSF1440-WW20 by6@carbonate.uits.iu.edu:/N/slate/by6/GSF1440_1441
scp GSF1441-WW20 by6@carbonate.uits.iu.edu:/N/slate/by6/GSF1440_1441


scp GSF1440-WW205-2_S5_R1_001_trimmed.fq.gz by6@carbonate.uits.iu.edu:/N/slate/by6/GSF1440_1441
#GSF1440-WW205-3_S9_R1_001_trimmed.fq.gz
#GSF1440-WW206-2_S6_R1_001_trimmed.fq.gz
#GSF1440-WW206-3_S10_R1_001_trimmed.fq.gz

#GSF1441-WW205-1_S5_R1_001_trimmed.fq.gz
#GSF1441-WW206-1_S6_R1_001_trimmed.fq.gz

GSF1440-WW205

# Unzip all files pasted 

gunzip GSF2848*	


################################
## Generate STAR Genome Index ##
################################
	
# Make a directory to store the genome files in the GSF2848 folder 
	
mkdir -p genome


# Then go back to /N/slate/by6/GSF2848 Always go back when you run any alignment etc. 

	
# Download and unpack the genome assembly.
	
curl $ASSEMBLY | gunzip > ./genome/assembly.fasta
	

# Download and unpack the genome annotation.
	
curl $ANNOTATION | gunzip > ./genome/annotation.gtf
	

# Create a directory to store the index.
	
mkdir -p genome/index
	

# Create the STAR genome index.
# --genomeSAindexNbases 12 was recommended by software.
	
	STAR \
	  --runThreadN 4 \
	  --runMode genomeGenerate \
	  --genomeDir genome/index \
	  --genomeFastaFiles genome/assembly.fasta \
	  --sjdbGTFfile genome/annotation.gtf \
	  --genomeSAindexNbases 12
	

###########################
## Align Reads to Genome ##
###########################
	
# Create an output directory for aligned reads.
	
mkdir -p results/aligned
	
  
# Align the reads.
	
FASTQ=$GSF2848*
	
  
	for FASTQ in ${FASTQ[@]}; do
	  PREFIX=results/aligned/$(basename $FASTQ .fastq)_
	  STAR \
	    --runThreadN 8 \
	    --outFilterMultimapNmax 1 \
	    --outFilterScoreMinOverLread .66 \
	    --outFilterMismatchNmax 10 \
	    --outFilterMismatchNoverLmax .3 \
	    --runMode alignReads \
	    --genomeDir genome/index \
	    --readFilesIn $FASTQ \
	    --outFileNamePrefix $PREFIX \
	    --outSAMattributes All \
	    --outSAMtype BAM SortedByCoordinate
	done

##Following has to be done additionally_20210817

# Indexing the BAM files.

BAMS=($(find ./results/aligned -name "*\.bam"))

for BAM in ${BAMS[@]}; do
  samtools index $BAM
  done
  
#remove genome* and results* files
  
  ####################
## Count Features ##
####################

# Create an output directory for read counts.

mkdir -p results/counts

# Count reads.

BAMS=$(find ./results/aligned -name "*\.bam")

featureCounts \
  -a genome/annotation.gtf \
  -o results/counts/counts.tsv \
  -t gene \
  -g gene_id \
  --largestOverlap \
  --readExtension3 150 \
  --primary \
  -s 2 \
  -T 8 \
  ${BAMS}
  
done
