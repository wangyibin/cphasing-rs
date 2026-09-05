# Library of [C-Phasing](https://github.com/wangyibin/CPhasing.git)

[![Bioconda Downloads](https://anaconda.org/bioconda/cphasing-rs/badges/downloads.svg)](https://anaconda.org/channels/bioconda/packages/cphasing-rs/overview)

## Installation

- Install with Conda ([Bioconda package](https://anaconda.org/channels/bioconda/packages/cphasing-rs/overview))
```
conda install bioconda::cphasing-rs
```



- Download from release
```
wget https://github.com/wangyibin/cphasing-rs/releases/download/latest/cphasing-x86_64-unknown-linux-musl.tar.gz
tar xzvf cphasing-x86_64-unknown-linux-musl.tar.gz 

```

- From source code
```
git clone https://github.com/wangyibin/cphasing-rs.git

cd cphasing-rs
pixi run install
```

## Optimize orientation

`cphasing-rs optimize` defaults to `--orientation-method banded-legacy`, restoring
the historical banded orientation and signed-block refinement pipeline. It uses a
rank window of 3, sqrt-links weights, at least 3 links per pair, and an input-sign
prior of 0.05. Confidence filtering and physical-span limits are disabled. After
orientation, it considers reverse-complement blocks of up to 32 contigs for at
most 4 accepted moves, with a relative local-band gain threshold of 0.0001. Each
accepted block is followed by another orientation pass.

Use `--orientation-block-span 0` to disable blocks, or
`--orientation-block-passes 2` to reduce the refinement budget. To select the
previous conservative policy, use `--orientation-method banded`; its defaults
remain a 0.95 confidence threshold, a 0.05 flip-bp fraction, and disabled blocks.
`banded-legacy` rejects confidence or bp limits because its internal orientation
passes do not implement those gates. `--orientation-method legacy` selects the
older ALLHiC orientation algorithm, which is distinct from `banded-legacy`.

The restored defaults reproduce the 2026-09-04 n500k alfalfa Hi-C result; this
does not establish optimal settings for other datasets. The experimental
`robust` and `banded-contact` modes remain opt-in.
