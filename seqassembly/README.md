# De Bruijn Graph Sequence Assembly
This module implements a partial assembly of a DNA sequence given a input of reads.
Specifically, it will obtain the contigs and contig lengths from a given set
of reads by building a de Bruijn graph.

The program conducts a simple error detection that filters out any reads that
contains K-mers that occur only once in all reads.

The module can either be imported as a package to use the
DeBruijnGraph class and class methods or run directly from the command
line.

## Command Line Usage
The general format of running the module from the command line is as follows:

```sh
python3 assembler.py [args] < [reads_input_file]
```

### Arguments
* [args] can be any of the following arguments:

| Argument             | Usage                 | About                                                                                  |
| -------------------- | --------------------- | -------------------------------------------------------------------------------------- |
| Help                 | -h, --help            | show help message and exit                                                             |
| K-mer length         | --k                   | The length of each k-mer, default 31                                                   |
| Minimum contig length| --min_length          | Minimum length desired of contigs to be returned, default 0                            |
| Output directory     | --out_dir             | The desired directory to create output files, by default the current working directory.|

The following command can be run for more information on arguments:
```sh
python3 assembler.py -h
```

Help:
```sh
$ python3 align.py -h
usage: assembler.py [-h] [--k K] [--min_length MIN_LENGTH] [--output_dir OUTPUT_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --k K                 The length of each k-mer
  --min_length MIN_LENGTH
                        Minimum length desired of contigs to be returned
  --output_dir OUTPUT_DIR
                        The desired directory to create output files, by default 
                        the current workingdirectory.
```

### Input file

* [reads_input_file] is the file path of the input file used to import the
squences reads to assemble the DNA sequence.

## Examples
### Example 1
```sh
python3 assembler.py < tests/input/sequence_reads
```
Uses the default arguments to determine the good reads, contigs, and contig 
lengths. Outputs each as a text file to the current working directory

### Example 2
```sh
python3 assembler.py --k 31 --min_length 100 --output_dir ./tests/output < tests/input/sequence_reads
```
Explicitly specifies K-mer length, although it is the same as the default and
the desired minimum contig length which is now 100. Also specifies a directory
for output files.

## Testing
The unit testing framework can be run with the following terminal command from
any parent directory of `assembler/tests`:
```sh
pytest
```
