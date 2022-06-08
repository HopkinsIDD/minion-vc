## Vibrio cholerae typing from Oxford Nanopore data

This is a guide for taking reads generated on a Oxford Nanopore sequencer and mapping them against specific genes (_ctxA_, _wbeO1_, _wbfO139_, _tcpA_ classical, _tcpA_ El Tor) for the purpose of determining sample serotype. Note that this guide is not optimized for reads generated with short amplicons.

## Setting up your working directory

This script assumpes the following:

* That you have already set up your working directory as described in the [assembly instructions](https://github.com/HopkinsIDD/minion-vc/blob/master/assembly_instructions.md) and that it contains a `demux` directory containing filtered FASTQ files (`*.filt.fastq.gz`).

* That the `scripts` folder inside your working directory contains the `map_genes.sh` script.

* That all five gene-specific reference FASTAs are contained in a directory called `typing`. If you cloned this git repository, these reference files can be found in `ref_genomes/typing`.

Once you have the required files in your working directory, navigate to this directory and create a new folder where you will store your typing results:

```
cd my-project
```

```
mkdir typing_results
```

## Running the typing script

In order to map reads from your samples to the reference genes identified above, you will first need to create a file that contains the names of the samples you would like to run the typing script on. For example, this file should look something like this:

```
sample1
sample2
sample6
sample7
```

Save this file as `samplenames.txt` inside of your `my-project/typing_results` directory. This file can include all samples, even those (especially those) that did not produce complete assemblies.

Then run the typing script as follows:

```
bash scripts/map_genes.sh typing_results/samplenames.txt my-project path-to-typing-directory min_cov typing_results/
```

Where `my-project` is the absolute path to your project directory, and `path-to-typing-directory` is the path to the directory that contains the `typing` folder containing the gene-specific FASTA files (usually `my-project/ref_genomes`). Additionally `min_cov` is the minimum number of reads that must map to a specific position to consider that site "covered". We suggest using 20 as a good starting value.

Note: this script may take 10-20 minutes to run per sample.


## Evaluating typing results

The typing script will produce one output file (`samplename_typing_summary.txt`) per sample included in `samplenames.txt`. This file will look something like this:

```
gene    reads_mapped    mean_depth      frac_covered
ctxA    192     178.609 1.0000000
tcpA_classical  237     215.929 1.0000000
tcpA_eltor      234     220.079 1.0000000
toxR    187     171.393 .9921436
wbeO1   1270    259.135 1.0000000
wbfO139 819     60.042  .2467412
```

These files may be easier to view in a spreadsheet software such as Excel.

Genes found in this sample will have a fraction covered (`frac_covered`) of near 1.00, indicating that reads map to the entire gene. Genes not found in this sample will have a lower fraction covered. In the example above, we can say that this sample contains the _ctxA_, _wbeO1_, _tcpA_ classical, and _tcpA_ El Tor genes, but not the _wbfO139_ gene. We can therefore feel confident saying that this is O1 _Vibrio cholerae_, though we cannot definitively say if the genome belongs to the Classical or El Tor biotype.
