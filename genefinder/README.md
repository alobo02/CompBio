# Content Statistics Gene Finding
This module uses content statistics to help determine the approximate start and 
stop location of the SARS-CoV-2 spike protein codon sequence in the reference
viral genome.

This is done by obtaining the sequences of the spike protein and 2000bp windows
in the reference viral genome. Then, for each pair of spike and window sequences,
a difference distance metric is computed by summing the difference of codon usage
frequencies squared between the spike sequence and the window sequence over
all possible (64) codons.

These difference distance metrics are then plotted to observe how the distances
change at different window start positions in the viral genome.

## Command Line Usage
The general format of running the module from the command line is as follows:

```sh
python3 covid_content_stats.py [args]
```

### Arguments
* [args] can be any of the following arguments:

| Argument             | Usage                 | About                                                                                  |
| -------------------- | --------------------- | -------------------------------------------------------------------------------------- |
| Help                 | -h, --help            | Show help message and exit                                                             |
| Spike sequence fasta | --spike_file          | Path to the input spike sequence fasta file                                            |
| Viral Genome fasta   | --viral_file          | Path to the input viral genome fasta file                                              |
| Output directory     | --out_dir             | The desired directory to create output files, by default the current working directory |

The following command can be run for more information on arguments:
```sh
python3 covid_content_stats.py -h
```

Help:
```sh
$ python3 covid_content_stats.py -h
usage: covid_content_stats.py [-h] [--spike_file SPIKE_FILE] [--viral_file VIRAL_FILE] [--output_dir OUTPUT_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --spike_file SPIKE_FILE
                        input file for the spike protein sequence
  --viral_file VIRAL_FILE
                        input file for the viral reference genome
  --output_dir OUTPUT_DIR
                        The desired directory to create output files, by default the current working directory.
```

## Examples
### Example 1
```sh
python3 covid_content_stats.py --spike_file tests/input/spike.fasta --viral_file tests/input/fullgenome.fasta
```
Outputs distance text file and plot file to the current working directory.
