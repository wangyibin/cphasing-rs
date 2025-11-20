## [v0.2.2] - 20251120
- Switch cargo version from 1.81 to 1.91, which increased the performance of software

## [v0.2.1] - 20251119
- Used in CPhasing v0.2.7.beta.r304

## [v0.2.0] - 20250327
- First release to the public, which is used in CPhasing v0.2.6.r297


## [v0.0.16]
## Enhancement
- `aligner`, change to banded_semi_global_SW algorithm 


## [v0.0.15]
## New features
- `realign`, first method of realign

## [v0.0.14]
## Enhancement
- `contacts`, advance normalization method
- `kprune`, remove self contig pairs

## [v0.0.13]
## New features
- `pairs2clm`, convert pairs to clm file
- `count_re`, count the restrict enzyme site counts in genome

## [v0.0.12]
## Enhancement
- `kprune`, introduce parallel 

## [v0.0.11]
## New features
- `pairs2bam`, convert 4DN pairs to a pseudo bam 
- `kprune`,  a function to identify the allelic and cross-allelic contig pairs

## [v0.0.10]
## Enancement
- `simulator`, add porec simulator which modified reads by vcf

## [v0.0.9]
## Bug fixed
- `aligner`, fixed bug of the orientation of mod reads

## [v0.0.8]
## Enhancement
- `simulation`, simulation data from long read
## Bug fixed
- `modbam2fq`, fixed the bug of read orientation

## [v0.0.7]
## Enhancement
- `aligner`, use command to execute minimap2

## [v0.0.6]
## New features
- `aligner`, realign alignments by mod(memory leak)

## [v0.0.5]
## Enhancement 
- `modbam2fq`,  `modfa`, enhanced

## [v0.0.4]
### New features
- `modbam2fq`, convert modbam to fastq with moddified bases
- `modfa`, modify fasta by methyation bed.
## [v0.0.3]
### Enhancement 
- `optimize`
    - fixed the orderidx error
    - advanced the optimize that can optimize 20 contigs
### Bug fixed
- `fastx`, `get_chrom_sizes` newlines error
