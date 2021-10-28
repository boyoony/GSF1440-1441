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

scp GSF1440-WW20* by6@carbonate.uits.iu.edu:/N/slate/by6/GSF1440_1441
scp GSF1441-WW20* by6@carbonate.uits.iu.edu:/N/slate/by6/GSF1440_1441


# Unzip all files pasted 

gunzip GSF1440-WW205-2_S5_R1_001_trimmed.fq.gz
gunzip GSF1440-WW205-3_S9_R1_001_trimmed.fq.gz
gunzip GSF1440-WW206-2_S6_R1_001_trimmed.fq.gz
gunzip GSF1440-WW206-3_S10_R1_001_trimmed.fq.gz

gunzip GSF1441-WW205-1_S5_R1_001_trimmed.fq.gz
gunzip GSF1441-WW206-1_S6_R1_001_trimmed.fq.gz


################################
## Generate STAR Genome Index ##
################################
	
# Make a directory to store the genome files in the GSF1440_1441 folder 
	
mkdir -p genome


# Then go back to /N/slate/by6/GSF1440_1441 Always go back when you run any alignment etc. 

	
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
	
FASTQ=$GSF144*
	
  
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




############################
#Downloading SAILOR Program
############################

#In Carbonate
#Make directory for analysis

cd /geode2/home/u010/by6/Carbonate

mkdir GSF1440_1441
	

#Load the singularity module

module load singularity



############################
#Moving all files needed for SAILOR to Carbonate (not in a directory) 
#############################

#Go to genome directory 

cd /N/slate/by6/GSF1440_1441/genome


#copy and paste assembly.fasta files to Carbonate 

cp assembly.fasta /geode2/home/u010/by6/Carbonate


#Go to GSF1440_1441 folder

cd /N/slate/by6/GSF1440_1441


#Make a folder for SAILOR

mkdir sailor 


#Go to aligned directory

cd /N/slate/by6/GSF1440_1441/results/aligned 


# create merged bam file (only if you have multiple replicates) 

module load samtools


samtools merge -@ 8 ../../sailor/GSF1440-WW205-2_S5_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam GSF1440-WW205-3_S9_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam GSF1441-WW205-1_S5_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam

samtools merge -@ 8 ../../sailor/GSF1440-WW206-2_S6_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam GSF1440-WW206-3_S10_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam GSF1441-WW206-1_S6_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam



#Go to Carbonate 

cd /geode2/home/u010/by6/Carbonate


#Copy and paste merged bam files 

cp /N/slate/by6/GSF1440_1441/sailor/* .



#########################
#Get C. elegans reference genome file with SNPs
#########################

#On mac desktop, go to HundleyLab server - protocols - Bioinformatics - Bioinformatics notes - Files for Sailor 

#Copy and paste this to your desktop - c.elegans.WS275.snps.sorted.bed

#In mac desktop terminal, copy the file to Carbonate 

scp /Users/yang/Desktop/c.elegans.WS275.snps.sorted.bed by6@carbonate.uits.iu.edu: 
	

#In Carbonate (not in a directory), you should have bam files, snps.sorted.bed file, assembly.fasta file




############################
#Mark downloaded program as executable within directory
############################

chmod +x sailor-1.0.4
	
./sailor-1.0.4
	
y


#######
#Open yaml file with vim and edit filenames
#######

# open .yaml file
	
vim ce11_example.yaml
	
i
	


# update filenames using arrow keys
#GSF1440-WW205-2_S5_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam
#GSF1440-WW206-2_S6_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam
#assembly.fasta
#c.elegans.WS275.snps.sorted.bed


# You can put only one file each. If you have multiple files, perform SAILOR multiple times. 
	
# save and quit (ESC - shift+ZZ)
	

# rename .yaml file
	
mv ce11_example.yaml GSF1440-WW205.yaml
	

# run sailor
	
./sailor-1.0.4 GSF1440-WW205.yaml


#Repeat this for another file

cp GSF1440-WW205.yaml GSF1440-WW206.yaml


#####################
#Open yaml file with vim and edit filenames
#####################

# open .yaml file
	
vim GSF1440-WW206.yaml
	
i
	

# update filenames using arrow keys

#GSF1440-WW206-2_S6_R1_001_trimmed.fq_Aligned.sortedByCoord.out.bam
	
# save and quit (ESC - shift+ZZ)
	

# run sailor
	
./sailor-1.0.4 GSF1440-WW206.yaml



#####################
#Interpreting the SAILOR data 
#####################

#In the results folder, 

less -s GSF1440-WW205-2_S5_R1_001_trimmed.fq_Aligned.sortedByCoord.out.fwd.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed


#Secure copy and paste onto mac desktop
#In the mac desktop, 

scp by6@carbonate.uits.iu.edu:/geode2/home/u010/by6/Carbonate/GSF1440-WW205/results/*\.bed ~/Desktop/

scp by6@carbonate.uits.iu.edu:/geode2/home/u010/by6/Carbonate/GSF1440-WW206/results/*\.bed ~/Desktop/


#Open in excel file 


#Annotation 





scp /Users/yang/Desktop/GSF1440* by6@carbonate.uits.iu.edu:/N/slate/by6/GSF1440_1441
