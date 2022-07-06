## Phylogenetic analysis of O1 Vibrio cholerae

This is a guide for building a phylogenetic tree from _Vibrio cholerae_ O1 genomes. The steps below assume the following:

* That you have a set of _Vibrio cholerae_ whole genome sequences you would like to visualize on a phylogenetic tree (for instructions on how to assemble whole genome sequences from Oxford Nanopore data, see [these assembly instructions](https://github.com/HopkinsIDD/minion-vc/blob/master/assembly_instructions.md)).

* That your genomes were assembled to the N15961 _Vibrio cholerae_ O1 reference genome, accessions AE003852/AE003853 (available [here](https://github.com/HopkinsIDD/minion-vc/blob/master/ref_genomes/vc_reference.fasta)). Genomes assembled to any other reference will need to be aligned prior to phylogenetic analysis.

* That you have already masked known recombination regions in these genomes, for example using the gff file available [here](https://github.com/HopkinsIDD/minion-vc/blob/master/recomb_mask.gff). The masking process is included as part of the [assembly pipeline](https://github.com/HopkinsIDD/minion-vc/blob/master/assembly_instructions.md).

* That you have the software [`gubbins`](https://github.com/nickjcroucher/gubbins) installed on your computer (Croucher et al. 2015).

## Table of Contents

* [Selecting a background dataset](#selecting-a-background-dataset)
* [Concatenating sequences](#concatenating-sequences)
* [Recombination masking with Gubbins](#recombination-masking-with-gubbins)
* [Building a maximum likelihood tree](#visualizing-a-maximum-likelihood-tree)
* [Visualizing your phlogenetic tree](#visualizing-your-phylogenetic-tree)

## Selecting a background dataset

Before building a phylogenetic tree, it can be important to select a set of previously published sequences to compare your published genomes to. This can help you identify which lineages your sequences belong to, and to determine if newly generated sequences are similar to previously published ones.

A list of accession numbers covering the majority of publicly available _Vibrio cholerae_ O1 sequences is available [here](https://github.com/HopkinsIDD/minion-vc/blob/master/global_accessions_updated20220609.csv). These sequences can be dowloaded from [NCBI](https://www.ncbi.nlm.nih.gov/) in either raw data, contig, or assembled genome form. If downloading the genomes as raw data or contigs, they will need to be assembled to the same O1 reference genome (N15961) before proceeding. You will also need to ensure these sequences are masked using the same file/method as used for the new sequences you want to add to the tree.

Please reach out to swohl@scripps.edu to inquire about _Vibrio cholerae_ alignments and assembled genomes from publicly available data, or to suggest inclusion of additional sequences in the global accessions file above.

## Concatenating sequences

To complete this step, you will need:

* The path to a folder containing newly assembled (and masked) genomes
* The path to a file containing (masked) background genomes you want to add to your phylogenetic tree

Before starting, make sure you are in whatever directory you are using as your **working directory**:

```
cd my-project
```

Where `my-project` is the path to your project directory. You may also want to make a new directory within your project directory to store phylogenetic results:

```
mkdir trees
```

Once you have set up your working directory, concatenate all of these sequences into one file as follows: 

```
cat path-to-genomes-directory/*.mask.fasta path-to-background-dataset > trees/all_sequences.fasta
```

This should create a new file in the `trees` directory that contains all of the background sequences as well as all of the newly generated sequences you want to include in the tree. All of these sequences **must** have been assembled to the same _Vibrio cholerae_ O1 reference genome and masked using the same file and method.

## Recombination masking with Gubbins

You will now need to run a special recombination masking software called Gubbins. In addition to the sites manually masked previously, this software will automatically detect additional regions of the genome that show evidence of past recombination and therefore must be removed (or masked) prior to phylogenetic analysis.

Please note that running Gubbins on large bacterial datasets make require substantial (up to 64GB) memory and runtime (up to 24 hours). This step is often best performed on a compute cluster. You can skip this step for the purposes of making a rough tree that shows general relationships between sequences, but this step should be performed prior to publication and final interpretation of phylogenetic data.

Gubbins can be run as follows, where `num-threads` is the number of threads to use for processing, and `path-to-concatenated-fasta` is the path to the file containing all sequences to be included in the tree, created as described above:

```
cd my-project/trees

run_gubbins.py --verbose --threads num-threads --prefix gubbins path-to-concatenated-fasta
```

Running Gubbins will generate a number of output files. The file ending in `.filtered_polymorphic_sites.fasta` contains your sequences with recombination sites removed. Additionally, this file contains only polymorphic sites (sites that are the same across all genomes have been removed), as these are the sites necessary for building a maximum likelihood tree.

Sites removed by Gubbins due to evidence of recombination can be found in the file ending in `recombination_predictions.gff`.

## Building a maximum likelihood tree

IQ-TREE is a commonly used tool for building maximum likelihood trees and is freely available [here](http://www.iqtree.org/). It can be run from the command line or from the [IQ-TREE web server](http://iqtree.cibiv.univie.ac.at/), though the web server may struggle with datasets containing more than a few hundred sequences.

Before running IQ-TREE from the command line, make sure you are in whatever directory you are using as your **working directory**. If you have a specific directory for tree-building, navigate there as well:

```
cd my-project/trees
```

IQ-TREE can then be run from the command line as follows. This command may take several hours to run. If it takes more than 12 hours, you may have non-O1 sequences in your tree that should be removed before proceeding (follow the [typing instructions here](https://github.com/HopkinsIDD/minion-vc/blob/master/typing_instructions.md) to determine if sequences are O1 _Vibrio cholerae_).

```
iqtree -s gubbins.filtered_polymorphic_sites.fasta -o '5174_7_5' -nt AUTO -m TEST -bb 1000
```

This command will run IQ-TREE on the Gubbins output file, `gubbins.filtered_polymorphic_sites.fasta`. Update this file name if you would like to build a tree using a different file.

The `-o 5174_7_5` option tells IQ-TREE to root the tree on a sequence called A6 (accesion: ERR025382) originally published by Mutreja et al. Nature 2011. This sequence is a commonly used root for O1 _Vibrio cholerae_ trees containing sequences assembled to N15961, but the option can be omitted to generated an unrooted tree.

The `-nt AUTO` command asks the computer to specific the number of threads to use for computation, and `-m TEST` asks the software to test various nucleotide substitution models to use in the tree (this option may add 1-2 hours of runtime). You can also specify a specific nucleotide model to use with `-m GTR` or another similar model.

`-bb 1000` tells the IQ-TREE software to perform 1000 boostrap iterations.

IQ-TREE will generate several output files in the directory containing the input file you specified. The file ending in `.treefile` is the one that can be most easily used for visualization.

## Visualizing your phlogenetic tree

There are many ways to visualize phylogenetic trees. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) is one commonly used software for visualizing phylogenetic trees. You can open your newly generated tree in this software and annotate it to better understand the relationship between _Vibrio cholerae_ lineages in the tree.