## Vibrio cholerae Genome Assembly from Oxford Nanopore Data

This is a guide for taking reads generated on a Oxford Nanopore sequencer and assembly by aligning reads to the N16961 reference. Note that this guide is not optimized for reads generated with short amplicons.

## Table of Contents

* [Setting up your working directory](#setting-up-your-working-directory)
* [Installation](#installation)
* [Basecalling](#basecalling)
* [Visualizing run metrics](#visualizing-run-metrics)
* [Demultiplexing](#demultiplexing)
* [Assembly](#assembly)
* [Masking](#masking)

## Setting up your working directory

Before you start, you will need to clone this repository. Cloning the repository will ensure you have all the necessary scripts and folders to run the analysis.

First, let's create a project directory to store all the analyses for this sequencing run. This will be your **working directory** for the rest of this process. Open up a terminal window and type the following command:

```
mkdir my-project
```
Where **my-project** can be any name you chose to represent this sequencing run.

To enter this directory, type:

```
cd my-project
```

To get the full path to this directory, now type:

```
pwd
```

The output, which might look something like `/home/username/my-project/`, will be your **working directory** for this analysis. Anytime you are asked to provide the path to your working directory, copy in this information.

From inside your **working directory**, type the following commands to copy this github repository onto your computer:

```
apt update
apt install git
git clone https://github.com/HopkinsIDD/minion-vc.git .
```

#### Finding your raw data

Some of the commands below require you to know where your raw MinION data is stored. When you started your MinION run, you selected a data folder in which to store the results. The raw data is in a folder called `fast5` within this folder. The path to your raw data probably looks something like this:

```
name-of-data-folder/run-name/library-name/flowcell-name/fast5
```

Figure out the path to your raw data and provide it below whenever the command calls for the **path to raw data**.


#### Determining the number of processors on your computer

It can be helpful to know how many processors are available on your computer, since using more processors can make some scripts run faster. Run this command to figure out how many you have available:

```
cat /proc/cpuinfo | grep processor
```
Use the bottom number on your screen as the number you provide any time a script askes you for **num-cores**.


## Installation

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

If you have not installed the software needed for assembly, run the following installation script:

```
sudo bash scripts/install-software.sh path-to-software
```

If some of the software is already installed on your computer, make sure that the `path-to-software` matches where the software is currently installed. If this is the first time you are installing assembly software, we reccomend the following path:

```
~/software
```
Please note that you will need an internet connection to complete this step, and that downloading the software will take some time (approximately 30-60min, depending on the speed of your internet connection). If the script is stopped for any reason, you can safely restart it with the command above and it will continue where it left off. If you are ever asked to respond to a yes/no questions, type "y" to continue.


## Basecalling

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

The first step in analyzing MinION data is turning the raw signal into nucleotides (e.g., A,C,T,G). We call this process *basecalling* and we use a software called `guppy` to do it. This process can take several hours to complete, depending on how much data you have generated.

```
guppy_basecaller -r -v -c dna_r9.4.1_450bps_fast.cfg -i path-to-raw-data -s basecalled -x auto
```

## Visualizing run metrics

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

After basecalling, you can optionally visualize some key metrics from your run like this. This may take a few minutes, depending on how much data you generated:

```
mkdir nanoplot
NanoPlot --summary basecalled/sequencing_summary.txt -o nanoplot

```
After this command finishes, open a file browser and find the nanoplot folder inside your **working directory**. Double click on the file ending in `.html` to see a visualization of your sequencing run.


## Demultiplexing

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

Demultiplexing is the process of splitting up the reads in your sequencing run by sample. In the lab, you assigned each sample a unique barcode. We will use these barcodes to assign each read to a sample, essentially dividing up all of the MinION reads into the samples they originally came from.

Open a text editor on your computer (e.g., Notepad). We will type up the barcodes so they can be used in our code. Type up your barcodes in the format below. The white space between the sample name and the barcode number must be a tab character, and you should replace **NB##** with the appropriate barcode number (e.g., NB02) and **sample1-name**, etc. with the names of your samples.

```
sample1-name	NB##
sample2-name	NB##
sample3-name	NB##
```

(repeat for all samples on your MinION run)

Then save this file inside your **working directory** and name it `samplesheet.txt`.

Now go back to your command line window. If you look inside your project directory, you should see the sample sheet file you just created:

```
ls my-project
```

Now you can start the demultiplex. Keep in mind that this process can take several hours to days, depending on the speed of your computer and the amount of data you generated. If the process ever gets interrupted, you can safely run the command below again and it will pick up where it left off:

```
bash scripts/batch_demux.sh my-project samplesheet.txt num-cores
```

To determine what number to enter instead of **num-cores**, see [Determining the number of processors on your computer](#determining-the-number-of-processors-on-your-computer) above.

## Assembly

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

We are now going to assemble the reads from the MinION sequencer into (hopefully) a Vibrio choelrae genome. To do this, we use a software called “nanopolish” that improves our basecalls and identifies variants from the reference genome. You will need to run this software once for each sample. Run the software as below, replacing **sample-name** with the name of the sample you want to run:

```
bash scripts/call_variants.sh sample-name my-project path-to-raw-data num-cores
```

To determine the path to enter instead of **path-to-raw-data**, see [Finding your raw data](#finding-your-raw-data) above. To determine what number to enter instead of **num-cores**, see [Determining the number of processors on your computer](#determining-the-number-of-processors-on-your-computer) above.

When this finishes, create your final genome by running the following command:

```
python3 scripts/make_consensus.py sample-name my-project
```
This script also prints out the mean and median coverage of your genome and creates a coverage plot. You can open this coverage plot by opening a file browser and navigating to your **working directory**. Inside the `genomes` folder, there should be a file called `sample-name.coverage.pdf`. Double click this file to view the plot.

After running the above two commands on all your samples, their final genomes can be found in the `genomes` folder inside your **working directory**.

## Masking

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

In this step, we will use a previously-published list of potential recombinant sites to mask areas of the genome that could mess up our phylogenetic analysis. The `recomb_mask.gff ` file used to mask sites was originally published at https://figshare.com/s/d6c1c6f02eac0c9c871e by Weill et al.

Before running the masking script, we need to make a text file containing the full paths to all genomes we want to mask. Use the following command to create this file:

```
ls -d genomes/*.ref.fasta > genomes/genome_dirs.txt
```

Now we can run the masking script with the following command:

```
bash scripts/mask_fasta.sh my-project genomes/genome_dirs.txt recomb_mask.gff "mask" "gff_seqname"
```

These genomes now have recombinant regions masked, and can be used in phylogenetic analyses. For best results, run gubbins on the final, masked alignment prior to phylogenetic analysis.
