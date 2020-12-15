# repTime
RepTime calculates the replication timing of WGS data.

RepTime combined the function of read depth calculation, GC correction, and replication timing smooth together in single-step command.

## Install
```
git clone https://github.com/deepomicslab/RepTime RepTime
```
## Requirements
```
SamTools
Python3

Python packages:
NumPy
Pandas
sklearn
csaps
```
## Usage
```
usage: RepTime.py [-h] [-i INPUT] [-c CHR] [-r REF] [-u UNIQUE]
                  [-w WINDOW] [-o OUTDIR] [-p PROCESS] [-s STEP]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input bam file
  -c CHR, --chr CHR     chr for analysis.For example, chr1 or
                        chr1-3,chr5,chr14 or 1-4,7,8-10,22
  -r REF, --ref REF     faidx indexed reference sequence file for input bam
  -u UNIQUE, --unique UNIQUE
                        uniquely mapped region
  -w WINDOW, --window WINDOW
                        the size of sliding window, step size, and cut-off
                        weight in this window. Default value is 10k window
                        with step size of 2k, 90 percent region of this windon
                        was sequenced
  -o OUTDIR, --outdir OUTDIR
                        the directory of output result
  -p PROCESS, --process PROCESS
                        the number of process
  -s STEP, --step STEP  the step of repTime.1:mpileup; 2:GC; 3:sliding window;
                        4:filter; 5:smooth
```

## Example
```
python RepTime.py -i your_input.bam -r hg19.fasta
```
