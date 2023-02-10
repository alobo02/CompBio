# NeedlemanWunch
This module uses the Needleman-Wunch algorithm to determine a global DNA 
sequence alignment using a simple similarity scoring 
scheme with a linear gap penalty: 
* M for each matched pair in the alignment
* m for each mismatch
* g for each gap

The module can either be imported as a package to use the
GlobalSequenceAlignment class and class methods or run directly from the command
line.

## Command Line Usage
The general format of running the module from the command line is as follows:

$ python3 align.py [args] < [FASTA_input_file] > [output_file]

### Arguments
* [args] can be any of the following arguments:
  -h, --help            show this help message and exit
  -M MATCH_SCORE, --match_score MATCH_SCORE
                        Match score to be used for global alignment, default 4
  -m MISMATCH_SCORE, --mismatch_score MISMATCH_SCORE
                        Mismatch score to be used for global alignment, 
                        default -2
  -g GAP_PENALTY, --gap_penalty GAP_PENALTY
                        Linear gap penalty to be used for global alignment,
                        default, -2

Run the following command for more information on arguments:
$ python3 align.py -h

### Input file

* [FASTA_input_file] is the file path of the input file used to import the
squences to be aligned

### Outputfile

* [output_file] is the desired file path of the file to print the optimal global
aligments to

## Examples
$ python3 align.py < tests/input/aligntest.input1 > tests/output/my.output1
$ python3 align.py -M 3 -m -2 -g -2 < tests/input/aligntest.input2 > tests/output/my.output2
$ python3 align.py --match_score 4 --gap_penalty -1 < tests/input/aligntest.input3 > tests/output/my.output3