# Readme file for fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py
# Copyright (c) 2015 Thaddeus D. Seher (@tdseher). All rights reserved.
# Barcode, spacer, and primer sequences from Fadrosh et al. Microbiome 2014

# Use the following command to view the program usage:
python "fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py" -h

# Use this command to debarcode the sample data:
python "fadrosh-16s-simple-debarcoder-with-mismatches-config-single+paired.py" -p parameters.txt -i sample-r1.fastq sample-r2.fastq -o debarcoded-r1.fastq debarcoded-r2.fastq -b bsp-r1.fastq bsp-r2.fastq -f fail-r1.fastq fail-r2.fastq

# If using your own data, be sure to modify the "parameters" file to include the
# list of forward and reverse barcodes, spacers, and primers you need to debarcode.
# Additionally, you will need a text file containing sample names and their (idealized) 
# barcode sequences.

# After debarcoding, I suggest you use either the default Illumina 
# demultiplexing program, or give FASTQ-MULTX from the EA-UTILS package a try:
# (https://github.com/ExpressionAnalysis/ea-utils)

# FASTQ-MULTX program usage can be found here:
# (https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMultx.md)

# You can use this code to demultiplex the sample data:
fastq-multx -H -B barcodes.txt debarcoded-r1.fastq debarcoded-r2.fastq -o %-r1.fastq %-r2.fastq
