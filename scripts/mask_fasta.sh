#! /bin/bash
set -e 
set -o pipefail

# Mask sequences with a given masking file

#### INPUT FILES ####

DIR=$1
FASTAS=$2 # file containing list of genomes to mask
REGIONS=$3 # regions to mask
PREFIX=$4 # name to include in filename
SEQNAME=$5 # chromosome name


### BEGIN MASK FASTA SCRIPT ###

while read line; do

	echo "new line: $line"
	HEADER=$(awk 'NR==1{print $1}' $line)
	SAMPLENAME="${HEADER:1}"
	echo $SAMPLENAME

	if [ ! -f $DIR/genomes/$SAMPLENAME.$PREFIX.fasta ]; then

		echo "creating masked file for $SAMPLENAME"

		# change the fasta header to match the chrom name in the GFF file
		var=">$SEQNAME"
		sed -i "1s/.*/$var/" $line

		# mask regions in fasta file
		bedtools maskfasta -fi $line -bed $REGIONS -fo $DIR/genomes/$SAMPLENAME.$PREFIX.fasta

		# change back the fasta header
		# and update fasta header in new file
		sed -i "1s/.*/$HEADER/" $line
		sed -i "1s/.*/$HEADER/" $DIR/genomes/$SAMPLENAME.$PREFIX.fasta
	
	fi

done < $FASTAS
