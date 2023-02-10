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

```sh
python3 align.py [args] < [FASTA_input_file] > [output_file]
```

### Arguments
* [args] can be any of the following arguments:

| Argument           | Usage                 | About                                                           |
| ------------------ | --------------------- | --------------------------------------------------------------- |
| Help               | -h, --help            | show help message and exit                                      |
| Match score (M)    | -M, --match_score     | Match score to be used for global alignment, default 4          |
| Mismatch score (m) | -m, --mismatch_score  | Mismatch score to be used for global alignment, default -2      |
| Gap penalty (g)    | -g, --gap_penalty     | Linear gap penalty to be used for global alignment, default, -2 |

The following command can be run for more information on arguments:
```sh
python3 align.py -h
```

Help:
```sh
$ python3 align.py -h
usage: align.py [-h] [-M MATCH_SCORE] [-m MISMATCH_SCORE] [-g GAP_PENALTY]

optional arguments:
  -h, --help            show this help message and exit
  -M MATCH_SCORE, --match_score MATCH_SCORE
                        Match score to be used for global alignment
  -m MISMATCH_SCORE, --mismatch_score MISMATCH_SCORE
                        Mismatch score to be used for global alignment
  -g GAP_PENALTY, --gap_penalty GAP_PENALTY
                        Linear gap penalty to be used for global alignment
```

### Input file

* [FASTA_input_file] is the file path of the input file used to import the
squences to be aligned

### Output file

* [output_file] is the desired file path of the file to print the optimal global
aligments to

Specifying the output file is optional. If it is not specified, then the sequence
alignments will be printed to the terminal.

## Examples
### Example 1
```sh
python3 align.py < tests/input/aligntest.input1 > tests/output/my.output1
```
The following will be printed to the terminal:
```sh
python3 align.py < tests/input/aligntest.input1
TC-CAAATAGAC
TCGCAAATATAC
```

### Example 2
```sh
python3 align.py < tests/input/aligntest.input1 > tests/output/my.output1
```
Same as example 1, except the output is printed to tests/output/my.output1

### Example 3
```sh
python3 align.py -M 3 -m -2 -g -2 < tests/input/aligntest.input2 > tests/output/my.output2
```
Specifies a custom scoring scheme that is different from all default values

### Example 4
```sh
python3 align.py --match_score 4 --gap_penalty -1 < tests/input/aligntest.input3 > tests/output/my.output3
```
Specifies a custom scoring scheme where the match and gap penalty is different
from the default, but the use of the default mismatch score is implied.