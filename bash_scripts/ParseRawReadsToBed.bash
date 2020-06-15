#!/bin/bash
# This script takes an sra file as input and runs the file through the SRA toolkit, trimmomatic, 
# bowtie, samtools, and bedtools to create a bed file.
# The first argument should be the input sra data, the second input should be the path to the fasta file of sequences for trimmomatic,
#   and the third argument should be the path to the basename for the bowtie2 index files.

# Get trimmomatic's jar path.
trimmomaticPath=$(dpkg -L trimmomatic | grep .jar$ | head -1)

# Get the name of the sra file from the user input.
dataName=${1%.*}
dataDirectory=${1%/*}
echo
echo "Working with $1"

# Create the names of all the intermediate and output files.
rawFastq="$dataName.sr.fastq.gz"
trimmedFastq="${dataName}_trimmed.fastq.gz"
bowtieSAMOutput="$dataName.sam.gz"
BAMOutput="$dataName.bam.gz"
finalBedOutput="$dataName.bed"

# Convert the data to SRA format
echo "Converting to fastq format..."
fastq-dump --gzip -O $dataDirectory $1

# Trim the data
echo "Trimming adaptors..."
java -jar $trimmomaticPath SE $rawFastq $trimmedFastq "ILLUMINACLIP:$2:2:30:10"

# Align the reads to the genome.
echo "Aligning reads with bowtie"
bowtie2 -x $3 -U $trimmedFastq -S $bowtieSAMOutput

# Convert from sam to bam.
echo "Converting from sam to bam..."
samtools view -b -o $BAMOutput $bowtieSAMOutput

# Convert to final bed output.
echo "Converting to bed..."
bedtools bamtobed -i $BAMOutput > $finalBedOutput