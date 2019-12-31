# SCYN: a local optimal single cell CNV profiling method using dynamic programming

SCYN: a local optimal single cell CNV profiling method using dynamic programming

## Pre-requirements
* python3
* numpy>=1.16.1
* pandas>=0.23.4,<0.24
* tasklogger>=0.4.0
* scipy>=1.3.0
* [SCOPE](https://github.com/rujinwang/SCOPE)


### install requirements
```Bash
pip install -r requirements.txt
```
To install R package SCOPE, please refer to the README of [SCOPE](https://github.com/rujinwang/SCOPE). SCYN integrates the SCOPE to get the cell-by-bin reads depth matrix and perform the normalization. SCYN mainly focuses on finding the optimal CNV segmentation profiling using dynamic programming.

## Installation

### Installation with pip
To install with pip, run the following from a terminal:
```Bash
pip install scyn
```

### Installation from Github
To clone the repository and install manually, run the following from a terminal:
```Bash
git clone https://github.com/xikanfeng2/SCYN.git
cd SCYN
python setup.py install
```

## Usage

### Quick start
The following code runs SCYN.

```Python
import scyn

# create SCYN object
scyn_operator = scyn.SCYN()

# call cnv
# bam_dir is the input bam directory and output_dir is the output directory
scyn_operator.call(bam_dir, output_dir)

# store cnv matrix to a csv file
scyn_operator.cnv.to_csv('your file name')
```

### Parameters
```Python
SCYN(seq='single-end', bin_len=500, ref='hg19', reg='*.bam', mapq=40, K=10, verbose=1)
```
Parameters

* seq : string, optional, default: single-end
    The reads type: single-end or paired-end

* bin_len : int, optional, default: 500
    The bin length, default is 500K

* ref : string, optional, default: hg19
    The reference genome version: hg19 or hg38

* reg : string, optional, default: *.bam
    The regular expression to match all BAM files in your input directory.
    For example, "*.bam" will match all BAM files ended with '.bam'

* mapq : int, optional, default: 40
    The mapping quality cutoff when calculating the reads coverage

* K : int, optional, default: 10
    The predifined change points number for all chromosomes


* verbose : `int` or `boolean`, optional, default: 1

    If `True` or `> 0`, print status messages

## Cite us

## Help
If you have any questions or require assistance using SCYN, please contact us with xikanfeng2@gmail.com.