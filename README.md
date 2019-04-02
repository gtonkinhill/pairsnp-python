# pairsnp-python

[![Travis-CI Build Status](https://travis-ci.com/gtonkinhill/pairsnp-python.svg?branch=master)](https://travis-ci.com/gtonkinhill/pairsnp-python)

## Installation

This is a Python3 version of pairsnp written with the help of Simon Frost

```
git clone https://github.com/gtonkinhill/pairsnp-python
cd pairsnp-python
pip3 install .
```

## Quick Start

The python version can be run from the python interpreter as

```
from pairsnp import calculate_snp_matrix, calculate_distance_matrix

sparse_matrix, consensus, seq_names = calculate_snp_matrix(fasta.file.name)
d = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)
```

alternatively if installed using pip it can be used at the command line as


```
pairsnp -f /path/to/msa.fasta -o /path/to/output.csv
```

additional options include

```
Program to calculate pairwise SNP distance and similarity matrices.

optional arguments:
  -h, --help            show this help message and exit
  -t {sim,dist}, --type {sim,dist}
                        either sim (similarity) or dist (distance) (default).
  -n, --inc_n           flag to indicate differences to gaps should be
                        counted.
  -f FILENAME, --file FILENAME
                        location of a multiple sequence alignment. Currently
                        only DNA alignments are supported.
  -z, --zipped          Alignment is gzipped.
  -c, --csv             Output csv-delimited table (default tsv).
  -o OUTPUT, --out OUTPUT
                        location of output file.
```
