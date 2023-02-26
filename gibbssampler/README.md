# Gibbs Sampler
This module implements a simplified variant of the Gibbs sampling algorithm to find a 
common sequence motif.

Unlike the typical algorithm, this variant chooses the best position in s*
(i.e. the window given a position with the highest score).

The module can either be imported as a package to use the
SimpleGibbsSampler class and class methods or run directly from the command
line.

## Command Line Usage
The general format of running the module from the command line is as follows:

```sh
python3 gibbs.py [args] < [FASTA_input_file] > [output_file]
```

### Arguments
* [args] can be any of the following arguments:

| Argument           | Usage                 | About                                                                          |
| ------------------ | --------------------- | ------------------------------------------------------------------------------ |
| Help               | -h, --help            | show help message and exit                                                     |
| Motif length       | --motif_len           | The fixed motif length, default 6                                              |
| Alphabet (Language)| --lang                | List characterizing the residue language symbols, default ['A', 'C', 'G', 'T'] |
| Background prob    | --lang_prob           | Background probability of each residue, typically 1/len(lang) , default 0.25   |
| Seed               | --seed                | Integer seed for random number generator, default, 20                          |

The following command can be run for more information on arguments:
```sh
python3 gibbs.py -h
```

Help:
```sh
$ python3 align.py -h
usage: gibbs.py [-h] [--motif_len MOTIF_LEN] [--lang LANG [LANG ...]]
                [--lang_prob LANG_PROB] [--seed SEED]

optional arguments:
  -h, --help            show this help message and exit
  --motif_len MOTIF_LEN
                        The fixed motif length
  --lang LANG [LANG ...]
                        List characterizing the residue language symbols
  --lang_prob LANG_PROB
                        Background probability of each residue, typically 1/len(lang)
  --seed SEED           Integer seed for random number generator
```

### Input file

* [FASTA_input_file] is the file path of the input file used to import the
squences to be used to determine the MSA of the best motif.

### Output file

* [output_file] is the desired file path of the file to print the iteration log to.

Specifying the output file is optional. If it is not specified, then the logs 
will be printed to the terminal.

## Examples
### Example 1
```sh
python3 gibbs.py < tests/input/Gibbs.short.fasta > tests/output/output.seed20
```
Uses the default arguments to determine the MSA of the best motif of the 
sequences stored in Gibbs.short.fasta

### Example 2
```sh
python3 gibbs.py --motif_len 6 --seed 14 < tests/input/Gibbs.fasta > tests/output/output.seed14
```
Explicitly specifies motif length, although it is the same as the default and
changes to integer used to seed the random number generator to 14.

## Testing
The unit testing framework can be run with the following terminal command from
any parent directory of `gibbssampler/tests`:
```sh
pytest
```
