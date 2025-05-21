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

### slapquant

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

#### Modus operandi

`slapquant` runs an aligner (BWA-MEM2 by default) to align the given RNA-Seq reads to the genome. The alignments are filtered for entries that are clipped at the start and/or the end of the read. The spliced leader sequence, if present, needs to be found immediately upstream of the start of the matched fragment. Therefore, a SLAS can be identified at the location of the first matching nucleotide immediately preceded by the proper spliced leader sequence. Similarly, the polyadenylation motif (at least six A bases) needs to immediately follow the end of the matched fragment, hence a PAS site is identified by the coordinate of the end of the matched fragment succeeded by the polyA motif.

Splice leader sequence detection (activated by default if no sequence is specified on the command line) works on the assumption that in an appropriately trimmed read file, the second most abundant clipped sequence at the start of a read will be the spliced leader sequence.

### slapassign

```
slapassign gene_models.gff slas_pas.gff > gene_models_with_slas_pas.gff
```
This takes the spliced leader acceptor sites and the polyadenylation sites in `slas_pas.gff` (the output of `slapquant` above) and assigns them to the gene models given in `gene_models.gff`. This obiously won't make sense if the gene models are based off a different genome than the original run of `slapquant` and `slapassign` will try to warn in this case.

The output of the `slapassign` run is GFF containing the original gene models and the SLAS / PAS sites assigned to them (by specifying a "Parent" attribute to the sites).

#### Modus operandi

In the case of SLAS, `slapassign` looks for coding sequence features (CDS) downstream of the SLAS (taking the strand into account). If a PAS site is closer than the next CDS or there is no CDS left on the sequence region, the SLAS site remains unassigned. For PAS, the algorithm looks for CDS upstream. Again, if a SLAS is closer, the PAS remains unassigned. 

### slapidentify

```
slapidentify genome.fasta reads1.fastq reads2.fastq reads3.fastq > sl_sequence.txt
```
This will run `slapidentify` with `genome.fasta` as the reference genome, using `reads1.fastq`, `reads2.fastq` and `reads3.fastq` as the input RNASeq reads. This uses the same spliced leader sequence identification strategy as running `slapquant` without the `-S` argument, but only returns the identified spliced leader sequence.

This is only intended for advanced users. Spliced leader identification by this method is highly sensitive to sequencing data quality, particularly trimming of adapter sequences. If you have any doubt, then use the existing known sequence for a species.

## API

`slapquant` can be used as a python module, however the internal API is likely still subject to change and documentation is still lacking.