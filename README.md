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
conservative policy, use `--orientation-method banded`; its defaults
remain a 0.95 confidence threshold, a 0.05 flip-bp fraction, and disabled blocks.
The confidence margin excludes the queried contig's own input-sign penalty,
so the 0.05 prior no longer makes the 0.95 threshold unreachable. After filtering,
the remaining changes are checked together against the input's regularized
score; a worse combination is rolled back and reported as `rejected joint changes`.
`banded-contact` remains a compatibility alias for this corrected solver.
`banded-legacy` rejects confidence or bp limits because its internal orientation
passes do not implement those gates. `--orientation-method legacy` selects the
older ALLHiC orientation algorithm, which is distinct from `banded-legacy`.

`--orientation-block-context-weight W` optionally adds long-range contacts to
block selection, with `W` between 0 and 1 (default 0, disabled). It requires
`banded-legacy` with blocks enabled and is experimental. At the start of block
refinement it freezes all complete contig-pair quartets beyond
`--orientation-window` with sufficient links.
The context score is the negative weighted sum of log distances between contig
midpoints, using `--orientation-pair-weight`. Candidates must improve the local
score and the combined score `local + W * context`, and are ranked by the
combined gain. Context alone may decrease if the local gain compensates it.
Moving a pair inside the local window does not remove it from the context. Midpoint scoring
is independent of signs, so orientation DP cannot undo this check. With no
qualifying context pairs, or with weight 0, historical selection is used exactly
unless joint candidate comparison is enabled below. Logs report the context score,
pair count, and number of rejected candidates. Agreement
with this contact score does not guarantee a lower true assembly error rate;
the check is disabled by default. The n500k alfalfa Hi-C comparison worsened
EditDistance at weights 0.001, 0.01, and 0.1. Keep weight 0 for the restored
behavior; positive weights are retained for experiments, not recommended as
a replacement for the validated defaults.

`--orientation-block-candidates K` enables experimental joint order/orientation
comparison (`K` from 1 to 16, default 0 for historical selection). It requires
`banded-legacy` with blocks enabled. Each pass retains the top K blocks passing
the existing preliminary local/combined gain gates, independently reorients each
from the same pass input, and selects the full post-orientation band score plus
the optional weighted context contribution. The winning tour, including its
computed signs, is applied directly. Both final local and combined gains must
pass the same threshold. The input-sign prior guides each candidate's orientation
DP as before; the ranking score itself excludes that penalty. Final-score ties
prefer shorter blocks, then earlier positions. Logs count reoriented candidates
and selections that differ from the preliminary winner. This uses up to K
orientation passes per sweep; it does not search candidates rejected by the
preliminary gates or optimize a reference-based error metric.
On frozen n500k alfalfa Hi-C inputs, K=8 reduced EditDistance from 13.2991%
to 13.2710%, but slightly reduced adjacency F1 and worsened the reference-order
control from 7.7417% to 7.7755%. K=2 was unchanged in EditDistance and K=4 was
worse. These results do not support changing the default from 0. A larger K
improves candidate coverage for one pass; different subsequent search paths
need not beat the historical final score within the same pass budget.

The restored defaults reproduce the 2026-09-04 n500k alfalfa Hi-C result; this
does not establish optimal settings for other datasets. The experimental
`robust` mode remains opt-in.

`--resume` requires an existing tour containing every input contig exactly once.
It reads the last nonempty line and retains unsigned names as forward-oriented
for compatibility. Missing, empty, incomplete, duplicate, or unknown-contig
tours are rejected before optimization. The original tour remains in place
until a complete result has been written successfully to a temporary file in
the same directory. Publication uses an atomic replacement, with the previous
contents saved to `.tour.sav`, then `.tour.sav.1`, `.tour.sav.2`, etc. Existing
backups are never overwritten. Failed input validation or optimization preserves
both the original tour and its backups.
