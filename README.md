# SCYN: Single cell CNV profiling method using dynamic programming

SCYN: Single cell CNV profiling method using dynamic programming


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

In command line:
```shell
usage: python run-scyn.py [-h] [options] -i input_bams_dir

SCYN: Single cell CNV profiling method using dynamic programming efficiently
and effectively

required arguments:
  -i, --indir   <str> the input bams directory (default: None)

optional arguments:
  -o, --outdir  <str> the output directory (default: ./)
  --seq           <str> the reads type: single-end or paired-end. (default:
                    single-end)
  --bin_len       <int> the bin length, default is 500K. (default: 500)
  --ref           <str> the reference genome version: hg19 or hg38.
                    (default: hg19)
  --reg           <str> the regular expression to match all BAM files in
                    your input directory. For example, ".bam" will match all
                    BAM files ended with '.bam'. (default: *.bam)
  --mapq          <int> the mapping quality cutoff when calculating the
                    reads coverage. (default: 40)
  --verbose       <int> If > 0, print log messages. (default: 1)
  -h, --help
```

In Python:
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

### SCYN attributes
```Python
scyn_operator = scyn.SCYN()
```
 - `scyn_operator.cnv` is the copy number variants matrix.
 - `scyn_operator.segments` is the segments for each chromosome.
 - `scyn_operator.meta_info` is the meta information of cells, include gini and ploidy.



### SCYN Output Format
The output of `SCYN` consits of two cnv files and one meta file. 

 - `cnv.csv`: with cell as row and bin as column. This file can be used as the input of Oviz-SingleCell CNV analysis.
 - `cnv_T.csv`: with bin as column and cell as row, it is the transpose matrix of `cnv.csv`. This file can be parse by popular R packages like [`ExpressionSet`](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) for downstream analysis.
 - `segments.csv` is the cnv segments information for each chromosome.
 - `meta.csv`: with cell as row, and meta information as column. The default meta information is:
   + `c_gini`: stores the gini coeficient of each cell.
   + `c_ploidy`: stores the mean ploidy of each cell, it is calculated from `cnv.csv` (not the one SCOPE provide).
   
   User can manually add extra cell meta information like 'cell_type', 'cluster', or 'group' for downstream analysis. Prefix `c` here denotes numeric continuous value. The absence of prefix `c` denotes category meta information like 'group' or 'cluster'.

### Parameters
```Python
SCYN(seq='single-end', bin_len=500, ref='hg19', reg='*.bam', mapq=40, verbose=1)
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


* verbose : `int` or `boolean`, optional, default: 1

    If `True` or `> 0`, print status messages

## Cite us

## Help
If you have any questions or require assistance using SCYN, please contact us with xikanfeng2@gmail.com.