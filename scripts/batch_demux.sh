#!/bin/bash
set -e 
set -o pipefail

# Split a directory containing raw .fastq files into subdirectories
# each containing a prespecified number of files
# Porechop will be run on each sequentially to prevent memory issues
# due to loading the entire concatenated fastq.gz file

#### INPUT FILES ####

DIR=$1 # project directory
SAMPLE_SHEET=$2 # sample sheet
NUM_CORES=$3 # number of processors available

#### VALUES THAT CAN BE CONFIGURED ####

FASTQ="$DIR/basecalled" # assumed path to raw fastq files
NUM_FILES=50 # number of fastq files per subdirectory


#### BEGIN DEMUX SCRIPT ####

date

# fist check that the provided working directory is a real directory
# with a subdirectory as defined above
if [ ! -d $DIR ]; then
	echo "ERROR: invalid directory provided"
	exit 1
fi

if [ ! -d $FASTQ ]; then
	echo "ERROR: no directory $FASTQ; please check directory structure"
	exit 1
fi

DIR=${DIR%/} # remove trailing slash from directory name if necessary

# make demux directory if it does not already exist
if [ ! -d $DIR/demux ]; then
	mkdir $DIR/demux
fi

#### SPLIT FILES INTO SUBDIRECTORIES ####

# if fastq files have not yet been sorted into bins
# also confirms the demux directory is empty
# since there will be no bins directory if the demux has already completed
if [ ! -d $FASTQ/bins/ ] && [ -z "$(ls -A $DIR/demux)" ]; then
		
	# make a new subdirectory
	mkdir $FASTQ/bins
	cd $FASTQ/bins/

	# loop through fastq files and sort them into appropriate bins
	COUNTER=0
	for f in $FASTQ/*.fastq; do
		d=bin_$(printf %03d $(($COUNTER/$NUM_FILES+1))) # determine the name of the directory
		mkdir -p $d # make the new directory
		mv $f $d # move the current fastq file to the directory
		COUNTER=$((COUNTER + 1)) # increase the counter
	done

	echo "moved all the files, beginning demux"

# if the fastq files have already been sorted (this is a restart of the demux script)
else
	echo "resuming demux"
fi

#### RUN PORECHOP ON EACH SUBDIRECTORY ####

# bin directories are removed after their demux is completed
# so we can just loop through the remaining directories

# first check if the bins directory exists
# if it does not exist (e.g. demux is complete), we can start concatenating files

if [ -d $FASTQ/bins/ ]; then

	for d in $FASTQ/bins/*/
	do
		b=${d::-1}
		b=${b##*/}
		echo "$b demux starting"
		porechop --format fastq.gz -b $DIR/demux/$b -i $FASTQ/bins/$b -t $NUM_CORES
		mv $FASTQ/bins/$b/* $FASTQ/ # add code to move files out of bins directory
		rmdir $FASTQ/bins/$b
	done

	# remove bins directory alltogether
	rmdir $FASTQ/bins/
fi

# at this point, the demux directory should be populated
# if it is not, suggest rerunning the demux script
# and remove the bins directory if empty to reduce confusion

if [ ! -z "$(ls -A $DIR/demux)" ]; then
	echo "demux complete, concatenating files"
else
	if [! -z "$(ls -A $FASTQ/bins)"]; then
		rmdir $FASTQ/bins
	fi
	echo "ERROR: error demuxing basecalled data; consider rerunning this script with updated arguments"
	exit 1
fi

#### CONCATENATE RESULTS FROM ALL BINS ####

# be in working directory
# just in case full path to sample sheet is not provided
cd $DIR

while read line; do
	b=$(echo $line | awk '{print $2}')
	echo "BARCODE $b"
	declare -a fq
	declare -a nobc
	
	# add barcode files in each directory to list for cat
	# if a barcode file does not exist in a particular directory, continue with warning
	for d in $DIR/demux/*/; do
		nobc+=( $d\none.fastq.gz ) # make sure to concatenate reads with no matching barcode
		if [ -f $d$b.fastq.gz ]; then
			fq+=( $d$b.fastq.gz )
		else
			echo "WARNING: $d does not contain barcode $b"
		fi
	done

	cat ${nobc[@]} > $DIR/demux/none.fastq.gz
	
	# check to see if fq array is empty
	# if it is, do not create a fastq file for that barcode
	# since no reads matched it
	if [ ! ${#fq[@]} -eq 0 ]; then
		cat ${fq[@]} > $DIR/demux/$b.fastq.gz
	fi

	unset fq
done < $SAMPLE_SHEET

echo "concatenating files complete, renaming files"

# rename output files according to sample sheet

while read line; do
	SAMPLENAME=$(echo $line | awk '{print $1}').fastq.gz
	BARCODE=$(echo $line | awk '{print $2}').fastq.gz
	
	# confirm that this barcode exists in the demux folder
	# if not, continue with warning
	if [ -f $DIR/demux/$BARCODE ]; then
		mv $DIR/demux/$BARCODE $DIR/demux/$SAMPLENAME
	else
		echo "WARNING: barcode $BARCODE was not found in basecalled reads"
	fi
done < $SAMPLE_SHEET

echo "batch demux completed successfully"
date