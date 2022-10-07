# Fadrosh 16s debarcoder #

Readme file for `fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py`
Copyright (c) 2015 Thaddeus D. Seher (@tdseher). All rights reserved.
Barcode, spacer, and primer sequences from [Fadrosh et al. Microbiome 2014](https://doi.org/10.1186/2049-2618-2-6).

Now updated for Python 3.

## Description ##
In order to overcome the early cycle fluorescence flooding in Illumina technologies when sequencing nearly-identical amplicons, "heterogeneity spacers" can be added to the read start termini of experimental libraries.

```text
        ... GCATAGCGA                          GGACGCATTGAAGCGGAACGCTG...
        ... TGTAGACCA T                        GGACGCATTGAAGCGGAACGCTG...
        ... AACCGTAAA GT                       GGACGCATTGAAGCGGAACGCTG...
        ... GGACTTTAG CGA                      GGACGCATTGAAGCGGAACGCTG...
        ... ATTATCCAT ATGA                     GGACGCATTGAAGCGGAACGCTG...
        ... TTGAGACCA TGCGA                    GGACGCATTGAAGCGGAACGCTG...
        ... AGTACGATT GAGTGG                   GGACGCATTGAAGCGGAACGCTG...
        ... CCAGTCAAT CCTGTGG                  GGACGCATTGAAGCGGAACGCTG...
[5' Linker] [Index 1] [Heterogeneity spacer 1] [Experimental sequence] [Heterogeneity spacer 2] [Index 2] [3' Linker]
```

If you do this, then you need a way of trimming out these "heterogeneity spacers" from your resultant FASTQ files. You can use this software to do so easily.

## Program usage ##
Use the following command to view the program usage:
```sh
python "fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py" -h
```
> ```text
> usage: fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py
>        [-h] -p PARAMETERS -i FASTQ [FASTQ ...] [-o FASTQ [FASTQ ...]] [-b FASTQ FASTQ]
>        [-f FASTQ [FASTQ ...]] [-t] [-a]
> 
> Python program to debarcode (but not demultiplex) Illumina sequences prepared through the
> custom 16S experiment. The program cycles through each sequence in pairs of reads
> simultaneously (or merged reads), first matching the observed barcode with the best candidate
> barcode for both reads (or both ends of the merged read), then matching the expected spacer
> and primer for both reads (or both ends of the merged read). If matches are found, then the
> sequences are printed to output files and the barcode, spacer, and primer are printed to bsp
> files. Otherwise they are printed to the fail files. Copyright (c) 2015 Thaddeus D. Seher
> (@tdseher). All rights reserved.
> 
> options:
>   -h, --help            show this help message and exit
>   -t, --true            Add true barcodes instead of idealized barcodes to sequence headers
>   -a, --autofield       If Illumina format headers detected, and the second field is missing,
>                         then add it.
> 
> required arguments:
>   -p PARAMETERS, --parameters PARAMETERS
>                         path to parameters file, specifying the read1 and read2 barcode,
>                         spacer, and primer sequences as well as the number of mismatches
>                         allowed in each
>   -i FASTQ [FASTQ ...], --input FASTQ [FASTQ ...]
>                         path to FASTQ read1 (and optionally read2) for reading input
>                         sequences
> 
> one or more arguments required:
>   -o FASTQ [FASTQ ...], --output FASTQ [FASTQ ...]
>                         path to FASTQ read1 (and optionally read2) for writing debarcoded
>                         sequences
>   -b FASTQ FASTQ, --bsp FASTQ FASTQ
>                         path to FASTQ read1 and read2 (left and right in merged) for writing
>                         barcode+spacer+primer sequences
>   -f FASTQ [FASTQ ...], --fail FASTQ [FASTQ ...]
>                         path to FASTQ read1 (and optionally read2) for writing non-debarcoded
>                         sequences
> ```

## Description of `parameters` file ##
Blank lines and lines starting with the `#` character are ignored.

Column description:
- The first (left-most) column is the read number (`1` or `2` permitted).
- The second column is barcode sequence
- The third column is the spacer sequence
- the fourth column is the primer sequence

For each read (left-most column), there can be a single line defining the number (if > 1) or the frequency (if < 1) of mismatches in that segment of each read.
Otherwise, there can be any number of lines describing barcodes, spacers, and primers to segregate.

## Example ##
Use this command to debarcode the sample data:
```sh
python "fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py" -p parameters.txt -i sample-r1.fastq sample-r2.fastq -o debarcoded-r1.fastq debarcoded-r2.fastq -b bsp-r1.fastq bsp-r2.fastq -f fail-r1.fastq fail-r2.fastq
```

## Notes ##
If using your own data, be sure to modify the `parameters.txt` file to include the list of forward and reverse barcodes, spacers, and primers you need to debarcode.
Additionally, you will need a text file containing sample names and their (idealized) barcode sequences.

After debarcoding, I suggest you use either the default Illumina demultiplexing program, or give FASTQ-MULTX ([documentation](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMultx.md)) from the [EA-UTILS](https://github.com/ExpressionAnalysis/ea-utils) package a try:

You can use this code to demultiplex the sample data:
```sh
fastq-multx -H -B barcodes.txt debarcoded-r1.fastq debarcoded-r2.fastq -o %-r1.fastq %-r2.fastq
```
