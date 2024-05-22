# slapquant

A simple tool to find spliced leader acceptor and polyadenylation sites in Discoba genomes from RNASeq data.

## Installation

`slapquant` requires
 * `python` version 3.11 or above
 * `bwa-mem` (or a compatible aligner, we recommend BWA-MEM2 due to its speed)
 * `gawk` (which should be installed by default on any mainstream Linux distribution)

 The easiest way to install `bwa-mem2` is via conda by running:
 ```
 conda install -c bioconda bwa-mem2
 ```
Please see the [Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on how to install conda.

Install `slapquant` by running:
```
pip install git+https://github.com/Wheeler-Lab/slapquant
```

## Usage

```
slapquant genome.fasta reads1.fastq reads2.fastq reads3.fastq > slas_pas.gff
```
This will run `slapquant` with `genome.fasta` as the reference genome, using `reads1.fastq`, `reads2.fastq` and `reads3.fastq` as the input RNASeq reads. `slapquant` is agnostic towards single- vs. paired-end reads and will treat every filename given after the genome file as a separate reads file. This command will then output a GFF file containing `SLAS` (spliced leader acceptor site) and `PAS` (polyadenylation site) entries:
```
...
###
chromosome1     slapquant       PAS     1869078 1869079 .       +       .       ID=chromosome1_PAS_16691;usage=1
###
chromosome1     slapquant       SLAS    1869391 1869392 .       +       .       ID=chromosome1_SLAS_5265;usage=1
###
chromosome1     slapquant       PAS     1869772 1869773 .       -       .       ID=chromosome1_PAS_84284;usage=1
###
...
```
Each entry has a `usage` attribute, which reflects how many times a read aligned to the given site.

## API

`slapquant` can be used as a python module, however the internal API is likely still subject to change and documentation is still lacking.