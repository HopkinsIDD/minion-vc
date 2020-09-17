#!/usr/bin/env python

""" Function for converting nanopolish VCFs to genomes:
    calls only variants of type "SNP"
    and imposes a minimum depth for calling a non-'N' base

    Adapted from: https://github.com/blab/zika-seq/blob/master/pipeline/scripts/margin_cons.py

    MUST BE RUN WITH PYTHON3

"""

import numpy as np
import pandas as pd

import sys
import vcf

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


#### COMMAND LINE ARGUMENTS ####

SAMPLENAME = sys.argv[1]
DIR = sys.argv[2] # directory with data for the sample

DIR = DIR.rstrip('/') # remove trailing slash if necessary

#### VALUES TO CHANGE PRIOR TO RUN ####

DEPTH_THRESHOLD = 20
REFERENCE=DIR+"/ref_genomes/Vibrio_cholerae_O1_biovar_eltor_str_N16961_v2.fasta" # path to reference genome

HQ_THRESHOLD=0.8 # threshold for a high quality genome

vcffile = DIR+"/nanopolish/"+SAMPLENAME+".nanopolish.filt.rename.vcf.gz"
depthfile = DIR+"/nanopolish/"+SAMPLENAME+".ref.depth.txt"
outdir = DIR+"/genomes/"+SAMPLENAME+".ref.fasta"
outplot = DIR+"/genomes/"+SAMPLENAME+".coverage.pdf"

#### BEGIN CONSENSUS GENOME ASSEMBLY ####

print("output file for sample")
print(SAMPLENAME)

# load the depth file
# need to adjust positions for genome not broken into chromosomes
depth = pd.read_csv(depthfile,delim_whitespace=True,header=None,names=["chrom","pos","depth"])
chr1 = depth[depth["chrom"]=="AE003852"]
chr2 = depth[depth["chrom"]=="AE003853"]

# store as dictionary to make indexing faster
chr1_keys = list(chr1["pos"])
chr1_values = list(chr1["depth"])
chr1 = dict(zip(chr1_keys,chr1_values))

chr2_keys = list(chr2["pos"])
chr2_values = list(chr2["depth"])
chr2 = dict(zip(chr2_keys,chr2_values))

# load the reference sequence
cons = ''

seq = list(SeqIO.parse(open(REFERENCE),"fasta"))[0]
cons = list(seq.seq.upper())
assert len(cons)==4033501

# iterate through positions in the reference sequence
# determine depth of the sample at this position

for n, c in enumerate(cons):

    pos = n+1 # adjust for 0-indexing

    # assume coverage is zero unless we can find the position in the depth file
    cov = 0

    # adjust for multiple chromosomes
    if pos <= 2961182:

        if pos in chr1.keys():
            cov = chr1[pos]

    else:

        pos = pos-2961182 # adjust to new chromosome

        if pos in chr2.keys():
            cov = chr2[pos]

    # determine if the coverage is above the threshold
    # change base to 'N' if below threshold
    if cov < DEPTH_THRESHOLD:
        cons[n] = 'N'

# save the modified reference genome to a temp file
seq = ''.join(cons)
new_record = SeqRecord(Seq(seq),id=SAMPLENAME,description="")
filepath = DIR + "/tmp/" + SAMPLENAME + ".modified.ref.fasta"
SeqIO.write(new_record, filepath, "fasta")


# open the vcf for this sample
vcf = vcf.Reader(filename=vcffile)

# initialize lists to store VCF information
# this is necessary to deal with duplicates in VCF
poslist = []
reflist=[]
altlist=[]

for record in vcf:

    if record.ALT[0] != '.':
        # variant call

        # ignore indels
        if len(record.REF)>1 or len(record.ALT[0])>1:
            continue

        # input VCF is already filtered on support values
        # so do not include code to filter on these values here

        CHROM=record.CHROM
        POS=record.POS
        ALT=record.ALT[0]
        REF=record.REF

        # adjust position for multiple chromosomes
        if CHROM=="AE003853":
            POS=POS+2961182

        # confirm that the reference allele matches the current consensus
        # skip position if the consensus has already been changed to 'N'
        if cons[POS-1] != str(REF):
            
            # this should only happen if consensus is 'N'
            # unless this position is a duplicate
            if cons[POS-1]=='N':
                continue
            
            else:
                assert (POS in poslist)
                
                # also confirm the REF and ALT alleles are correct
                idx = poslist.index(POS)
                assert reflist[idx]==REF
                assert altlist[idx]==ALT
                
                # if this is a true duplicate SNP, do nothing
                continue

        # after ruling out all other issues
        # assign alternate allele to the consensus genome
        cons[POS-1] = str(ALT)

        # save information to lists
        poslist.append(POS)
        reflist.append(REF)
        altlist.append(ALT)

m=0
for base in cons:
    if base=='N':
        m = m+1

# save high quality genomes
if float((4033501-m)/4033501)>0.8:

    print("high quality genome, coverage = ",float((4033501-m)/4033501))

    # save the genome to file
    seq = ''.join(cons)
    new_record = SeqRecord(Seq(seq),id=SAMPLENAME,description="")
    SeqIO.write(new_record, outdir, "fasta")

# save lower quality genomes with warning
else:
    print("low quality genome, coverage = ",float((4033501-m)/4033501))

    # save the genome to file
    seq = ''.join(cons)
    new_record = SeqRecord(Seq(seq),id=SAMPLENAME,description="")
    SeqIO.write(new_record, outdir, "fasta")


#### PLOT COVERAGE FOR THIS GENOME ####

# calculate and print the mean and median coverage based on samtools depth
med_depth = depth["depth"].median()
print("median depth of coverage:")
print(med_depth)

mean_depth = depth["depth"].mean()
print("mean depth of coverage:")
print(mean_depth)

# create plot
fig = plt.figure(figsize=(12,4),frameon=False)
ax1 = fig.add_subplot(111)

ax1.plot(depth["pos"],depth["depth"],color='lightgray')
#ax1.plot(depth["pos"],depth["depth"],linestyle="None",marker='o',markerfacecolor='white',markeredgecolor='gray')

ax1.set_xlim(0,len(depth.index))

ymax = depth["depth"].max()
ax1.set_ylim(0,ymax)

ax1.set_ylabel("Read depth",rotation="vertical")
ax1.set_xlabel("Position along genome")

ax1.yaxis.set_label_coords(-0.06,0.5)

#tick mark and axes settings
ax1.tick_params(axis='x',       #changes apply to the x-axis
            which='both',      #both major and minor ticks are affected
            bottom=True,       #ticks along the bottom edge are on
            top=False,         #ticks along the top edge are off
            labelbottom=True)  #labels along the bottom edge are on

ax1.tick_params(axis='y',       #changes apply to the y-axis
            which='both',      #both major and minor ticks are affected
            left=True,         #ticks along the left edge are on
            right=False,       #ticks along the right edge are off
            labelleft=True)    #labels along the left edge are on

ax1.xaxis.grid(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

plt.savefig(outplot,format="pdf")
