## Vibrio cholerae Genome Assembly from Oxford Nanopore Data

This is a guide for taking reads generated on a Oxford Nanopore sequencer and assembly by aligning reads to the N16961 reference. Note that this guide is not optimized for reads generated with short amplicons. Instructions for using this pipeline are available in `assembly_instructions.md` and `analysis_instructions.md`.

Please note this repository is no longer actively maintained.

## Features

This pipeline performs the following functions:

- [x] Basecalling (FAST5 to FASTQ) using [guppy](https://nanoporetech.com/software/other/guppy) (unsupported as of 2023)
- [x] Visualizing Oxford Nanopore run metrics using [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [x] Read demultiplexing using [porechop](https://github.com/rrwick/Porechop) (unsupported as of 2018)
- [x] Filtering of low quality reads using [Filtlong](https://github.com/rrwick/Filtlong)
- [x] Read alignment using [minimap2](https://github.com/lh3/minimap2)
- [x] Variant calling using [nanopolish](https://nanopolish.readthedocs.io/en/latest/index.html)
- [x] Recombination masking using [gubbins](https://github.com/nickjcroucher/gubbins)
- [x] Maximum-likelihood phylogenetic inference of processed samples and background dataset using [iqtree](https://github.com/iqtree/iqtree2)


## Installation

1. Clone the minion-vc repository:
```commandline
git clone https://github.com/HopkinsIDD/minion-vc.git
```

2. Move to the pipeline source directory:
```commandline
cd minion-vc/
```

3. Install software needed for assembly by running the following installation script:
```commandline
bash scripts/install-software.sh path-to-software
```

## Usage

1. Create a directory specifically for the batch of samples you would like to analyze (called a project directory).
```commandline
mkdir my-project
```

2. Before starting, make sure you are in your **project directory**:
```commandLine
cd my-project
```

3. Follow insturctions detailed in [`assembly_instructions.md`](https://github.com/HopkinsIDD/minion-vc/blob/f927c241090fa12183686906c596c98e70c4f948/assembly_instructions.md) to proceed with assembly of genomic data.