# Library of [C-Phasing](https://github.com/wangyibin/CPhasing.git)
 

## Installation

```
git clone https://github.com/wangyibin/cphasing-rs.git

cd cphasing-rs
cargo build --release 

```


## Command Usage
```bash
Phasing and scaffolding based on Pore-C or Hi-C data

Usage: cphasing-rs <COMMAND>

Commands:
  aligner          align methylation reads
  realign          rescue secondary alignment by high-quality alignments
  kprune           Identify the allelic and cross-allelic contig pairs
  splitbam         split bam by record number
  splitfastq       split fastq by record number
  simulator        simulating test data
  digest           digest genome by restriction enzyme, output restriction sites
  count_re         count restriction enzyme sites, only support single enzyme
  cutsite          cut restriction site on pore-c reads
  paf2porec        convert paf to pore-c table
  porec2pairs      convert pore-c table to pairs
  porec-merge      Merge multiple pore-c table file into single file
  porec-intersect  According a bed file to intersection a pore-c table.
  paf2pairs        convert paf to pairs
  pairs2contacts   calculate the contacts between contigs
  pairs2clm        convert pairs to clm file, which is used for `allhic optimize`
  pairs2mnd        convert pairs to mnd file
  pairs2bam        convert pairs to pseudo bam file
  pairs-intersect  According a bed file to intersection a pairs file.
  chromsizes       generate chromsizes file
  prunepairs       prune pairs
  modbam2fq        convert modified bam to fastq with modified base
  modfa            modify fasta by bedMethy file
  optimize         optimize contigs (developing)
  help             Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print versio
```
