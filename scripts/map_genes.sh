#! /bin/bash
set -e 
set -o pipefail

# Align reads for a given sample to the following genes:
## ctxA -> cholera toxin; associated with pandemic lineages O1 and O139
## wbeO1 -> O1 antigen
## wbfO139 -> O139 antigen
## tcpA classical
## tcpA El Tor


# Calculate metrics on mapped reads
## number and fraction of reads mapped
## depth of coverage
## percent coverage (positions with a minimum coverage depth)

date

# get sample name file and project directory from command line
NAMES=$1 # text file with sample names
DIR=$2 # path to project directory
REFDIR=$3 # path to reference genome directory
MIN_COV=$4 # minimum number of reads mapping to a site to count it as covered
OUTDIR=$5

# run minimap2 to to align reads from each sample to each typing reference
while read line; do

	# create a file to store the typing results from this sample
	if [ ! -f $OUTDIR/$line\_typing_summary.txt ]; then
		printf "gene\treads_mapped\tfrac_mapped\tmean_depth\tfrac_covered\n" > \
		$OUTDIR/$line\_typing_summary.txt
	fi

	echo $OUTDIR/$line\_typing_summary.txt

	# loop through all typing genes
	for ref in $REFDIR/typing/*; do

		# get the name of the gene we are testing
		name=${ref%%.fasta}
		name=${name##*/}

		echo $name

		# align reads with minimap2
		minimap2 -a -x map-ont $ref $DIR/demux/$line.filt.fastq.gz | \
		samtools sort - -o $OUTDIR/$line.$name.aln.sorted.bam

		# run samtools depth on the resulting file
		# so that we can calculate depth-related metrics
		samtools depth -a -m 10000 $OUTDIR/$line.$name.aln.sorted.bam > \
		$OUTDIR/$line.$name.depth.txt

		# calculate number of reads mapped
		reads_mapped=$(samtools view -c -F 4 $OUTDIR/$line.$name.aln.sorted.bam)

		# calculate fraction of reads mapped
		# first calculate number of reads
		reads=$(samtools view -c $OUTDIR/$line.$name.aln.sorted.bam)
		frac_mapped=$(echo "scale=7; $reads_mapped / $reads" | bc)

		# calculate mean depth
		mean_depth=$(awk '{ total += $3 } END { print total/NR }' $OUTDIR/$line.$name.depth.txt)

		# calculate percent of gene covered
		pos_covered=$(awk -v MIN_COV=$MIN_COV 'BEGIN {count = 0} $3 > MIN_COV {count++} END {print count}' $OUTDIR/$line.$name.depth.txt)
		total_pos=$(wc -l < $OUTDIR/$line.$name.depth.txt)
		frac_covered=$(echo "scale=7; $pos_covered / $total_pos" | bc)

		# put all data into a file
		printf "$name\t$reads_mapped\t$mean_depth\t$frac_covered\n" >> \
		$OUTDIR/$line\_typing_summary.txt

	done

done < $NAMES
