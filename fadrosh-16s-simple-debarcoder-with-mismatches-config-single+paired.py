#!/usr/bin/env python

# Copyright (c) 2015 Thaddeus D. Seher (@tdseher). All rights reserved.

# Define imports
import sys
import argparse
import string

# Define global variables
__author__ = "Thaddeus D. Seher"
__twitter__ = "@tdseher"
__date__ = "2015"

__description__ = """\
Python program to debarcode (but not demultiplex) Illumina sequences \
prepared through the custom 16S experiment. The program cycles through \
each sequence in pairs of reads simultaneously (or merged reads), first \
matching the observed barcode with the best candidate barcode for both reads \
(or both ends of the merged read), then matching the expected spacer and \
primer for both reads (or both ends of the merged read). If matches are \
found, then the sequences are printed to output files and the barcode, spacer, and primer \
are printed to bsp files. Otherwise they are printed to the fail files. \
Copyright (c) {__date__} {__author__} ({__twitter__}). All rights reserved. \
""".format(**locals())

# 1) search for barcode+spacer+primer sequence at left end of read
# 2) find the exact bp locations and split the read: barcode|spacer|primer|read
# 3) require strict matching in barcode+spacer
#    but allow the primer to have mismatches

def revcomp(dna):
    """Returns the reverse-complement of a string containing DNA characters"""
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]

def sequence_hamming_distance(first, second, max_differences):
    """Counts number of differences between two sequences. Both sequences
    are left-justified, and measured along the length of the first sequence
    only."""
    differences = 0
    for i in range(len(first)):
        if (first[i] != second[i]):
            differences += 1
    
    if (0 < max_differences < 1):
        if (differences > (len(first) * max_differences)):
            return None
        else:
            return differences
    else:
        if (differences > max_differences):
            return None
        else:
            return differences

def ambiguous_hamming_distance(first, second, max_differences):
    """Counts number of differences between two sequences. Both sequences
    are left-justified, and measured along the length of the first sequence
    only. IUPAC ambiguity characters may be used."""
    iupac = {
        #'a': ['a'],
        #'c': ['c'],
        #'g': ['g'],
        #'t': ['t'],
        #'r': ['a', 'g'],
        #'y': ['c', 't'],
        #'m': ['a', 'c'],
        #'k': ['g', 't'],
        #'w': ['a', 't'],
        #'s': ['c', 'g'],
        #'b': ['c', 'g', 't'],
        #'d': ['a', 'g', 't'],
        #'h': ['a', 'c', 't'],
        #'v': ['a', 'c', 'g'],
        #'n': ['a', 'c', 'g', 't'],
        
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'M': ['A', 'C'],
        'K': ['G', 'T'],
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
    }
    
    differences = 0
    for i in range(len(first)):
        m = [j for j in iupac[first[i]] if j in iupac[second[i]]]
        if (len(m) == 0):
            differences += 1
    
    if (0 < max_differences < 1):
        if (differences > (len(first) * max_differences)):
            return None
        else:
            return differences
    else:
        if (differences > max_differences):
            return None
        else:
            return differences

def assign_sequence(sequence, candidates, database, max_score=2):
    """Attempt to pair the experimental sequence's barcode with its best
    candidate barcode. Returns the index of the best candidate. If input
    sequence has multiple best candidates, then return None."""
    # candidates is a list of known barcodes that we expect
    # max_score is the maximum distance from a known barcode to assign
    # database holds a list of observed barcodes and their calculated scores
    
    match = None
    try:
        match = database[sequence]
    except KeyError:
        best_score = max_score
        best_candidates = []
        for i, c in enumerate(candidates):
            #score = align_sequence(sequence, c, best_score)
            score = sequence_hamming_distance(sequence, c, best_score)
            if (score != None):
                if (score == 0):
                    best_candidates = [i]
                    best_score = score
                    break
                elif (score < best_score):
                    best_score = score
                    best_candidates = [i]
                elif (score == best_score):
                    best_candidates.append(i)
        if (len(best_candidates) == 1):
            database[sequence] = best_candidates[0]
            match = best_candidates[0]
        else:
            database[sequence] = None
            match = None
        
    return match

def load_parameters(path):
    """Parse the file that specifies the read1 and read2 barcodes, spacers,
    and primers as well as their alignment models"""
    # If N == whole number, then this is the number of mismatches allowed:
    #  0 mismatches, 1 mismatches, 2 mismatches, etc
    # if 0 < N < 1, or a float, then this is the proportion if mismatches
    #  allowed. That is, if the sequence length is 10, and N = 0.2, then
    #  it may have at most 2 mismatches.
    parameters = {}
    barcodes = {}
    
    with open(path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                if not line.startswith("#"):
                    sline = line.split("\t")
                    try:
                        par_list = map(float, sline[1:])
                        parameters[sline[0]] =  par_list
                    except ValueError:
                        try:
                            barcodes[sline[0]].append(sline[1:])
                        except KeyError:
                            barcodes[sline[0]] = [sline[1:]]
    
    return parameters, barcodes

def check_pair(quad1, quad2):
    """Tests if the headers in two FASTQ sequences match"""
    return (quad1[0].split(" ", 1)[0] == quad2[0].split(" ", 1)[0])

def process_single(quad1, database1, database2, true=False, autofield=False):
    """Debarcode merged reads, concatenating their barcode sequences and
    adding them to the headers"""
    
    # we assume all barcodes are the same length
    seq1 = quad1[1]
    query1 = seq1[:len(barcodes["1"][0][0])]
    i1 = assign_sequence(query1, map(lambda x: x[0], barcodes["1"]), database1, parameters["1"][0])
    
    seq2 = revcomp(quad1[1])
    query2 = seq2[:len(barcodes["2"][0][0])]
    i2 = assign_sequence(query2, map(lambda x: x[0], barcodes["2"]), database2, parameters["2"][0])
    
    if ((i1 != None) and (i2 != None)):
        barcode1 = barcodes["1"][i1][0]
        barcode2 = barcodes["2"][i2][0]
        spacer1 = barcodes["1"][i1][1]
        spacer2 = barcodes["2"][i2][1]
        primer1 = barcodes["1"][i1][2]
        primer2 = barcodes["2"][i2][2]
        
        # Check the number of mismatches on both spacers and primers one at a time,
        # but breaking opon the first one to return None
        if (
            # make sure the sequence is long enough to cover all expected bases
            (len(seq1) >= (len(barcode1)+len(spacer1)+len(primer1) + len(barcode2)+len(spacer2)+len(primer2))) and
            # match the spacer sequences
            (sequence_hamming_distance(spacer1, seq1[len(barcode1):], parameters["1"][1]) != None) and
            (sequence_hamming_distance(spacer2, seq2[len(barcode2):], parameters["2"][1]) != None) and
            # match the primer sequences
            (ambiguous_hamming_distance(primer1, seq1[len(barcode1)+len(spacer1):], parameters["1"][2]) != None) and
            (ambiguous_hamming_distance(primer2, seq2[len(barcode2)+len(spacer2):], parameters["2"][2]) != None)
        ):
            
            new_quad1 = []
            
            if true:
                new_quad1.append(format_new_header(quad1[0], query1, query2, autofield))
            else:
                new_quad1.append(format_new_header(quad1[0], barcode1, barcode2, autofield))
            
            new_quad1.append(quad1[1][len(barcode1)+len(spacer1)+len(primer1):-1*(len(barcode2)+len(spacer2)+len(primer2))])
            new_quad1.append(quad1[2])
            new_quad1.append(quad1[3][len(barcode1)+len(spacer1)+len(primer1):-1*(len(barcode2)+len(spacer2)+len(primer2))])
            
            bsp_quad1 = []
            bsp_quad1.append(new_quad1[0])
            bsp_quad1.append(quad1[1][:len(barcode1)+len(spacer1)+len(primer1)])
            bsp_quad1.append(new_quad1[2])
            bsp_quad1.append(quad1[3][:len(barcode1)+len(spacer1)+len(primer1)])
            
            bsp_quad2 = []
            bsp_quad2.append(new_quad1[0])
            bsp_quad2.append(quad1[1][-1*(len(barcode2)+len(spacer2)+len(primer2)):])
            bsp_quad2.append(new_quad1[2])
            bsp_quad2.append(quad1[3][-1*(len(barcode2)+len(spacer2)+len(primer2)):])
            
            return [new_quad1], [bsp_quad1, bsp_quad2]
        
    return [None, None], [None, None]

def process_pair(quad1, quad2, database1, database2, true=False, autofield=False):
    """Debarcode the paired reads, concatenating their barcode sequences and
    adding them to the headers"""
    
    #sys.stderr.write("   quad1=" + str(quad1) + "\n")
    #sys.stderr.write("   quad2=" + str(quad2) + "\n")
    
    # we assume all barcodes are the same length
    query1 = quad1[1][:len(barcodes["1"][0][0])]
    i1 = assign_sequence(query1, map(lambda x: x[0], barcodes["1"]), database1, parameters["1"][0])
    
    query2 = quad2[1][:len(barcodes["2"][0][0])]
    i2 = assign_sequence(query2, map(lambda x: x[0], barcodes["2"]), database2, parameters["2"][0])
    
    #sys.stderr.write("      i1=" + str(i1) + "\n")
    #sys.stderr.write("  query1=" + query1 + "\n")
    #sys.stderr.write("      i2=" + str(i2) + "\n")
    #sys.stderr.write("  query2=" + query2 + "\n")
    
    if ((i1 != None) and (i2 != None)):
        barcode1 = barcodes["1"][i1][0]
        barcode2 = barcodes["2"][i2][0]
        spacer1 = barcodes["1"][i1][1]
        spacer2 = barcodes["2"][i2][1]
        primer1 = barcodes["1"][i1][2]
        primer2 = barcodes["2"][i2][2]
        
        #sys.stderr.write("barcode1=" + barcode1 + "\n")
        #sys.stderr.write("barcode2=" + barcode2 + "\n")
        #sys.stderr.write("\n")
        
        # Check the number of mismatches on both spacers and primers one at a time,
        # but breaking opon the first one to return None
        if (
            # make sure the sequence is long enough to cover all expected bases
            (len(quad1[1]) >= (len(barcode1)+len(spacer1)+len(primer1))) and
            (len(quad2[1]) >= (len(barcode2)+len(spacer2)+len(primer2))) and
            # match the spacer sequences
            (sequence_hamming_distance(spacer1, quad1[1][len(barcode1):], parameters["1"][1]) != None) and
            (sequence_hamming_distance(spacer2, quad2[1][len(barcode2):], parameters["2"][1]) != None) and
            # match the primer sequences
            (ambiguous_hamming_distance(primer1, quad1[1][len(barcode1)+len(spacer1):], parameters["1"][2]) != None) and
            (ambiguous_hamming_distance(primer2, quad2[1][len(barcode2)+len(spacer2):], parameters["2"][2]) != None)
        ):
            
            new_quad1 = []
            new_quad2 = []
            
            if true:
                new_quad1.append(format_new_header(quad1[0], query1, query2, autofield))
                new_quad2.append(format_new_header(quad2[0], query1, query2, autofield))
            else:
                new_quad1.append(format_new_header(quad1[0], barcode1, barcode2, autofield))
                new_quad2.append(format_new_header(quad2[0], barcode1, barcode2, autofield))
            
            new_quad1.append(quad1[1][len(barcode1)+len(spacer1)+len(primer1):])
            new_quad1.append(quad1[2])
            new_quad1.append(quad1[3][len(barcode1)+len(spacer1)+len(primer1):])
            
            new_quad2.append(quad2[1][len(barcode2)+len(spacer2)+len(primer2):])
            new_quad2.append(quad2[2])
            new_quad2.append(quad2[3][len(barcode2)+len(spacer2)+len(primer2):])
            
            bsp_quad1 = []
            bsp_quad1.append(new_quad1[0])
            bsp_quad1.append(quad1[1][:len(barcode1)+len(spacer1)+len(primer1)])
            bsp_quad1.append(new_quad1[2])
            bsp_quad1.append(quad1[3][:len(barcode1)+len(spacer1)+len(primer1)])
            
            bsp_quad2 = []
            bsp_quad2.append(new_quad2[0])
            bsp_quad2.append(quad2[1][:len(barcode2)+len(spacer2)+len(primer2)])
            bsp_quad2.append(new_quad2[2])
            bsp_quad2.append(quad2[3][:len(barcode2)+len(spacer2)+len(primer2)])
            
            return [new_quad1, new_quad2], [bsp_quad1, bsp_quad2]
        
    return [None, None], [None, None]

def format_new_header(header, barcode1, barcode2, autofield=False):
    split_header = header.split(" ")
    if (len(split_header) == 2):
        # Illumina format for second field in header
        #  Pos Element           Requirements  Description
        #  0   <read>            [12]          Read number. 1 can be single read or read 2 of paired-end
        #  1   <is filtered>     [YN]          Y if the read is filtered, N otherwise
        #  2   <control number>  [0-9]+        0 when none of the control bits are on, otherwise it is an even number.
        #  3   <index sequence>  [ACTG]*       Index sequence
        split_second_tag = split_header[1].split(":")
        if ((len(split_second_tag) == 4) and (split_second_tag[1] in ["Y", "N", "y", "n"])):
            split_second_tag[3] = barcode1 + barcode2
            return split_header[0] + " " + ":".join(split_second_tag)
    if autofield:
        # Illumina format for first field in header
        #  Pos Element           Requirements  Description
        #  0   @                 @             Each sequence identifier line starts with @
        #  0   <instrument>      [a-zA-Z0-9_]+ Instrument ID
        #  1   <run number>      [0-9]+        Run number on instrument
        #  2   <flowcell ID>     [a-zA-Z0-9]+
        #  3   <lane>            [0-9]+        Lane number
        #  4   <tile>            [0-9]+        Tile number
        #  5   <x_pos>           [0-9]+        X coordinate of cluster
        #  6   <y_pos>           [0-9]+        Y coordinate of cluster
        split_first_tag = split_header[0].split(":")
        if ((len(split_first_tag) == 7) and (split_first_tag[0][0] == '@')):
            split_header.insert(1, "1:N:0:" + barcode1 + barcode2)
            return " ".join(split_header)
    
    if (header[-1] == ";"):
        return header + "barcode=" + barcode1 + barcode2 + ";"
    elif (header[-1] == " "):
        return header + barcode1 + barcode2
    else:
        return header + " " + barcode1 + barcode2

def read_quad(flo):
    """Read four lines--a sequence entry, or "quad"--from an open FASTQ file"""
    quad = []
    for line in flo:
        quad.append(line.rstrip())
        if (len(quad) == 4):
            break
    return quad

# class argparse.Action
#  __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None, required=False, help=None, metavar=None)
# Action instances should be callable, so subclasses must override the __call__ method, which should accept four parameters:
#   parser - The ArgumentParser object which contains this action.
#   namespace - The Namespace object that will be returned by parse_args(). Most actions add an attribute to this object using setattr().
#   values - The associated command-line arguments, with any type conversions applied. Type conversions are specified with the type keyword argument to add_argument().
#   option_string - The option string that was used to invoke this action. The option_string argument is optional, and will be absent if the action is associated with a positional argument.
# The __call__ method may perform arbitrary actions, but will typically set attributes on the namespace based on dest and values.

class RequireSameAction(argparse.Action):
    """Action to require "-o/--output", "-b/--bsp", and "-f/--fail" arguments
    to have the same number of entries as the "-i/--input" argument"""
    def __call__(self, parser, namespace, values, option_string=None):
        if (len(values) > 2):
            raise argparse.ArgumentError(self, 'expected 1 or 2 arguments')
        
        if (len(values) != len(vars(namespace)['input'])):
            raise argparse.ArgumentError(self, 'expected same number of entries as argument -i/--input')
        
        #try:
        #    return super(RequireSameAction, self).__call__(parser, namespace, values, option_string)
        #except NotImplementedError:
        setattr(namespace, self.dest, values)

class RequirePairAction(argparse.Action):
    """Action to require either 1 or 2 entries for "-i/--input" argument"""
    def __call__(self, parser, namespace, values, option_string=None):
        if (len(values) not in [1, 2]):
            raise argparse.ArgumentError(self, 'expected 1 or 2 arguments')
        
        #try:
        #    return super(RequirePairAction, self).__call__(parser, namespace, values, option_string)
        #except NotImplementedError:
        setattr(namespace, self.dest, values)

def parse_arguments():
    """Creates the argument parser and stores positional parameters"""
    
    # Create the argument parser
    parser = argparse.ArgumentParser(
        description=__description__,
        epilog='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add mandatory arguments
    required_group = parser.add_argument_group("required arguments", description=None)
    required_group.add_argument("-p", "--parameters", required=True,
        help="path to parameters file, specifying the read1 and read2 barcode, \
             spacer, and primer sequences as well as the number of mismatches \
             allowed in each"
    )
    required_group.add_argument("-i", "--input",
        metavar="FASTQ", nargs="+", required=True,
        action=RequirePairAction,
        help="path to FASTQ read1 (and optionally read2) for reading input \
             sequences"
    )
    
    # Add partially-required arguments
    outputs_group = parser.add_argument_group("one or more arguments required", description=None)
    outputs_group.add_argument("-o", "--output",
        metavar="FASTQ", nargs="+", default=[],
        action=RequireSameAction,
        help="path to FASTQ read1 (and optionally read2) for writing \
             debarcoded sequences"
    )
    outputs_group.add_argument("-b", "--bsp",
        metavar="FASTQ", nargs=2, default=[],
        help="path to FASTQ read1 and read2 (left and right in merged) for \
             writing barcode+spacer+primer sequences"
    )
    outputs_group.add_argument("-f", "--fail",
        metavar="FASTQ", nargs="+", default=[],
        action=RequireSameAction,
        help="path to FASTQ read1 (and optionally read2) for writing \
             non-debarcoded sequences"
    )
    
    # Add optional arguments
    parser.add_argument("-t", "--true",
        action="store_true",
        help="Add true barcodes instead of idealized barcodes to sequence \
             headers"
    )
    parser.add_argument("-a", "--autofield",
        action="store_true",
        help="If Illumina format headers detected, and the second field is \
             missing, then add it."
    )
    
    # parse the arguments
    args = parser.parse_args()
    
    # Require at least one output type to be specified
    if (0 == len(args.output) == len(args.bsp) == len(args.fail)):
        parser.error('expected at least one -o/--output, -b/--bsp, or -f/--fail argument')
    
    return args

def main():
    """Run the program"""
    
    # load arguments and parse them
    args = parse_arguments()
    
    # load the parameters file
    global parameters, barcodes
    parameters, barcodes = load_parameters(args.parameters)
    
    # create dictionary lookup-tables to store barcode's best matches
    # each sequence read in will first be compared to the dictionary's keys
    # so that a full alignment can be skipped if it has already been hashed.
    database1 = {}
    database2 = {}
    
    # open input filehandles
    flo_in = []
    for f in args.input:
        flo_in.append(open(f, 'r'))
    
    # open output filehandles
    flo_output = []
    for f in args.output:
        flo_output.append(open(f, 'w'))
        
    # open bsp filehandles
    flo_bsp = []
    for f in args.bsp:
        flo_bsp.append(open(f, 'w'))
    
    # open fail filehandles
    flo_fail = []
    for f in args.fail:
        flo_fail.append(open(f, 'w'))
    
    # create list to store pairs of sequence lines
    quads = []
    for i in range(len(flo_in)):
        quads.append([])
    
    # Read in the first sequence from both input files
    for i in range(len(flo_in)):
        quads[i] = read_quad(flo_in[i])
    
    # cycle through the input FASTQ files simultaneously
    while (quads[0] != []):
        if (len(quads) == 1):
            new_quads, bsp_quads = process_single(quads[0], database1, database2, args.true, args.autofield)
        else:
            # make sure the sequence headers are properly paired, otherwise exit
            if not check_pair(quads[0], quads[1]):
                sys.stderr.write("headers not equal\n")
                sys.stderr.write("quad1=" + str(quads[0][0]) + "\n")
                sys.stderr.write("quad2=" + str(quads[1][0]) + "\n")
                sys.exit(1)
        
            # debarcode
            new_quads, bsp_quads = process_pair(quads[0], quads[1], database1, database2, args.true, args.autofield)
        
        # print debarcoded sequences
        if new_quads[0]:
            for i, flo in enumerate(flo_output):
                flo.write("\n".join(new_quads[i]) + "\n")
            
            for i, flo in enumerate(flo_bsp):
                flo.write("\n".join(bsp_quads[i]) + "\n")
        else:
            for i, flo in enumerate(flo_fail):
                flo.write("\n".join(quads[i]) + "\n")
        
        # read in the next sequence from both input files
        for i in range(len(flo_in)):
            quads[i] = read_quad(flo_in[i])

if (__name__ == '__main__'):
    main()
