#! /bin/bash
set -e 
set -o pipefail

#Use nanopolish to index reads based on raw signal data,
#then use minimap2 to align reads to a pre-specified reference,
#and finally call variants based on this alignment

#### INPUT FILES ####

# get sample name and directories from command line
NAME=$1 # sample name
DIR=$2 # path to project directory
READS=$3 # path to raw .fast5 files
NUM_CORES=$4 # number of cores available

#### CHECK DIRECTORIES ####

# fist check that the provided working directory is a real directory
# with reads for the desired sample in the demux directory
if [ ! -d $DIR ]; then
	echo "ERROR: invalid directory provided"
	exit 1
fi

if [ ! -f $DIR/demux/$NAME.fastq.gz ]; then
	echo "ERROR: reads for sample $NAME not in $DIR/demux/"
fi

echo "output file for $NAME"

DIR=${DIR%/} # remove trailing slash from directory name if necessary

# create necessary directories
# if they do not exist already
if [ ! -d $DIR/nanopolish ]; then
	mkdir $DIR/nanopolish
fi

if [ ! -d $DIR/tmp ]; then
	mkdir $DIR/tmp
fi

if [ ! -d $DIR/genomes ]; then
	mkdir $DIR/genomes
fi

#### VALUES TO CHANGE PRIOR TO RUN ####

READ_SUMMARY=$DIR/basecalled/sequencing_summary.txt # path to sequencing summary file generated by guppy

REF=$DIR/ref_genomes/vc_reference.fasta # path to reference genome
MAKERANGE=$DIR/scripts/nanopolish_makerange.py

REF_NAME_1="AE003852"
REF_NAME_2="AE003853"
REF_LEN_1=2961181
REF_LEN_2=1072318

SNP_SUPPORT=0.75
REQ_DEPTH=100

#### BEGIN NANOPOLISH PIPELINE ####

# filter input data if file sizes are large
date
echo "running filtlong"
if [ ! -f $DIR/demux/$NAME.filt.fastq.gz ]; then
        filtlong --keep_percent 90 --target_bases 800000000 $DIR/demux/$NAME.fastq.gz | gzip > $DIR/demux/$NAME.filt.fastq.gz
fi

#index reads individually for each sample
date
echo "running nanopolish index"
if [ ! -f $DIR/demux/$NAME.filt.fastq.gz.index ]; then
	nanopolish index -d $READS -s $READ_SUMMARY $DIR/demux/$NAME.filt.fastq.gz
fi

# map reads to reference
date
echo "running minimap2"
if [ ! -f $DIR/nanopolish/$NAME.aln.ref.sorted.bam ]; then
	minimap2 -ax map-ont -t $NUM_CORES $REF $DIR/demux/$NAME.filt.fastq.gz | \
		samtools sort -@ $NUM_CORES -o $DIR/nanopolish/$NAME.aln.ref.sorted.bam \
		-T $DIR/nanopolish/tmp/$NAME.reads.tmp
	samtools index $DIR/nanopolish/$NAME.aln.ref.sorted.bam
fi

echo "number of reads in sample $NAME:"
samtools view -c $DIR/nanopolish/$NAME.aln.ref.sorted.bam

echo "number of mapped reads in sample $NAME:"
samtools view -c -F 4 $DIR/nanopolish/$NAME.aln.ref.sorted.bam

date
echo "calling variants with nanopolish"

# use makerange function bundled with nanopolish to break up into window sizes
# this is to prevent loading all fast5 data into memory each time

# run makerange and save output to file
if [ ! -f $DIR/nanopolish/ref_ranges.txt ]; then
        python $MAKERANGE $REF > $DIR/nanopolish/ref_ranges.txt
fi

# loop through window sizes and run nanopolish on each
count=1

while read line; do
        if [ ! -f $DIR/nanopolish/$NAME.win$count.nanopolish.vcf ]; then
        
                #name output file using counter
                outfile="$DIR/nanopolish/$NAME.win$count.nanopolish.vcf"
                
                #run nanopolish here
                nanopolish variants -o $outfile \
                        -p 1 \
                        -w $line \
                        -q dam \
                        -r $DIR/00_demux/$NAME.filt.fastq.gz \
                        -b $DIR/nanopolish/$NAME.aln.ref.nanopolish.sorted.bam \
                        -g $REF \
                        -t $NUM_CORES
        fi

        echo "nanopolish for window $count complete"
        
        #increase counter used for naming files
        count=$(($count + 1))

done < $DIR/nanopolish/ref_ranges.txt

echo "nanopolish completed successfully for all windows"

#merge vcfs created by window size

#create file containing list of vcf files to merge
#also update the vcf header to contain info about the reference
#bcftools concat will not work otherwise
idx=1
touch $DIR/nanopolish/$NAME.vcf.files.txt
echo "##contig=<ID=$REF_NAME_1,length=$REF_LEN_1>" > $DIR/nanopolish/headertext.txt
echo "##contig=<ID=$REF_NAME_2,length=$REF_LEN_2>" >> $DIR/nanopolish/headertext.txt

while [ $idx -lt $count ]; do

        if [ ! -f $DIR/nanopolish/$NAME.win$idx.annotate.nanopolish.vcf ]; then

                #update vcf header
                bcftools annotate -h $DIR/nanopolish/headertext.txt -o $DIR/nanopolish/$NAME.win$idx.annotate.nanopolish.vcf \
                        --no-version $DIR/nanopolish/$NAME.win$idx.nanopolish.vcf

                #get outfile name and append to file
                echo "$DIR/nanopolish/$NAME.win$idx.annotate.nanopolish.vcf" >> $DIR/nanopolish/$NAME.vcf.files.txt
        fi

        #increase counter
        idx=$(($idx + 1))
done

if [ ! -f $DIR/nanopolish/$NAME.nanopolish.vcf ]; then
        #concatenate the vcf files for this sample
        bcftools concat --no-version -o $DIR/nanopolish/$NAME.nanopolish.vcf -f $DIR/nanopolish/$NAME.vcf.files.txt
fi

# filter VCF on hard-coded SNP support and read depth
if [ ! -f $DIR/nanopolish/$NAME.nanopolish.filt.vcf ]; then
	bcftools filter --no-version -i "INFO/SupportFraction>$SNP_SUPPORT" $DIR/nanopolish/$NAME.nanopolish.vcf | \
		bcftools filter --no-version -i "INFO/TotalReads>$REQ_DEPTH" -o $DIR/nanopolish/$NAME.nanopolish.filt.vcf
fi

#compress and index files for optional additional vcf manipulation
if [ ! -f $DIR/nanopolish/$NAME.nanopolish.filt.rename.vcf.gz ]; then
	
	bcftools view --no-version -O z -o $DIR/nanopolish/$NAME.nanopolish.filt.vcf.gz $DIR/nanopolish/$NAME.nanopolish.filt.vcf
	echo "sample $NAME" > $DIR/nanopolish/reheader-vcf.txt
	bcftools reheader -s $DIR/nanopolish/reheader-vcf.txt -o $DIR/nanopolish/$NAME.nanopolish.filt.rename.vcf.gz \
		$DIR/nanopolish/$NAME.nanopolish.filt.vcf.gz
	bcftools index --threads $NUM_CORES $DIR/nanopolish/$NAME.nanopolish.filt.rename.vcf.gz
fi

# calculate depth at each position
# needed to make consensus genome
samtools depth -a -m 10000 $DIR/nanopolish/$NAME.aln.ref.sorted.bam > $DIR/nanopolish/$NAME.ref.depth.txt

echo "nanopolish script completed successfully"