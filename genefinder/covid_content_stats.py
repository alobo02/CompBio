#!/usr/bin/env python3
"""Compute for codon usage statistics and gene finding for the canonical 
SARS-CoV-2 virus "surface glycoprotein" (i.e. "spike protein").

Command Line Usage
------------------
The general format of running the module from the command line is as follows:

$ python3 covid_content_stats.py [args]

where

[args] can be any of the following arguments:
  -h, --help            show this help message and exit
  --spike_file SPIKE_FILE
                        Input file for the spike protein sequence
  --viral_file VIRAL_FILE
                        Input file for the viral reference genome
  --output_dir OUTPUT_DIR
                        The desired directory to create output files, by default the current working directory.

Run the following command for information on arguments:
$ python3 covid_content_stats.py -h

Usage Examples
--------------
$ python3 covid_content_stats.py --spike_file tests/input/spike.fasta --viral_file tests/input/fullgenome.fasta

Author: Alexander A. Lobo
"""
from ast import arg
import typing
from typing import Union
import sys
import argparse
import numpy as np
import pathlib
from matplotlib import pyplot as plt

CODONS = {
    'TCA',    # Serina
    'TCC',    # Serina
    'TCG',    # Serina
    'TCT',    # Serina
    'TTC',    # Fenilalanina
    'TTT',    # Fenilalanina
    'TTA',    # Leucina
    'TTG',    # Leucina
    'TAC',    # Tirosina
    'TAT',    # Tirosina
    'TAA',    # Stop
    'TAG',    # Stop
    'TGC',    # Cisteina
    'TGT',    # Cisteina
    'TGA',    # Stop
    'TGG',    # Triptofano
    'CTA',    # Leucina
    'CTC',    # Leucina
    'CTG',    # Leucina
    'CTT',    # Leucina
    'CCA',    # Prolina
    'CCC',    # Prolina
    'CCG',    # Prolina
    'CCT',    # Prolina
    'CAC',    # Histidina
    'CAT',    # Histidina
    'CAA',    # Glutamina
    'CAG',    # Glutamina
    'CGA',    # Arginina
    'CGC',    # Arginina
    'CGG',    # Arginina
    'CGT',    # Arginina
    'ATA',    # Isoleucina
    'ATC',    # Isoleucina
    'ATT',    # Isoleucina
    'ATG',    # Methionina
    'ACA',    # Treonina
    'ACC',    # Treonina
    'ACG',    # Treonina
    'ACT',    # Treonina
    'AAC',    # Asparagina
    'AAT',    # Asparagina
    'AAA',    # Lisina
    'AAG',    # Lisina
    'AGC',    # Serina
    'AGT',    # Serina
    'AGA',    # Arginina
    'AGG',    # Arginina
    'GTA',    # Valina
    'GTC',    # Valina
    'GTG',    # Valina
    'GTT',    # Valina
    'GCA',    # Alanina
    'GCC',    # Alanina
    'GCG',    # Alanina
    'GCT',    # Alanina
    'GAC',    # Acido Aspartico
    'GAT',    # Acido Aspartico
    'GAA',    # Acido Glutamico
    'GAG',    # Acido Glutamico
    'GGA',    # Glicina
    'GGC',    # Glicina
    'GGG',    # Glicina
    'GGT',    # Glicina
}
# Source: https://gist.github.com/juanfal/09d7fb53bd367742127e17284b9c47bf


def count_codon(seq: str) -> dict:
    """Counts the occurrence of ALL possible codons in given sequence.

    Parameters
    ----------
    seq : str
        Input sequence

    Returns
    -------
    count_dict : dict
        Dictionary of counts for each possible codon
    """
    # Initialize dict with 0 count for each codon
    count_dict = {codon_key: 0 for codon_key in CODONS}

    # Iterate through valid codons in sequence and increment counts
    for codon in get_codons(seq):
        count_dict[codon] += 1

    return count_dict


def get_codons(seq: str) -> list:
    """Gets list of valid codons from a sequence."""
    return [seq[i:i+3] for i in range(0, len(seq), 3) if is_valid_codon(seq[i:i+3])]


def is_valid_codon(codon: str) -> bool:
    """Check if a given codon is valid."""
    return codon in CODONS


def compute_window_distances(genome_seq: str, spike_dict: dict) -> list[tuple[int, float]]:
    """Compute the window difference distances for each possible window of 
    length 2000 in the genome.

    Parameters
    ----------
    genome_seq : str
        Reference viral genome to window
    spike_dict : dict
        Dictionary of codon counts for spike protein sequence

    Returns
    -------
    list[tuple[int, float]]
        List of tuples in the form (window start position, window distance)
    """
    # Define window length
    window_length = 2000

    # Initialize lists to store window positions affiliated distances
    positions = []
    distances = []

    # Iterate through windows in viral genome and compute distance metric
    for pos in range(len(genome_seq)-window_length):
        # Get window associated condon_dict
        window = genome_seq[pos:pos+window_length]
        window_dict = count_codon(window)

        # Update lists
        positions.append(pos)
        distances.append(compute_distance(spike_dict, window_dict))

    return list(zip(positions, distances))


def compute_distance(spike_dict: dict, window_dict: dict) -> float:
    """Compute the difference distance metric between the spike protein sequence 
    and a given window from the full viral genome.

    Formula given by

    D = Σ_c (F_S - F_W)**2

    where,

    F_S = S_c / N_S
    N_S = Σ_c Sc
    F_W = W_c / N_W
    N_W = Σ_c W_c

    and
    S_c : Number of occurrences of codon c in the true spike sequence
    W_c : Number of occurrences of codon c in that window

    Parameters
    ----------
    spike_dict : dict
        Dictionary of codon counts for spike protein sequence
    window_dict : dict
        Dictionary of codon counts for given sequence window

    Returns
    -------
    distance : float
        Difference distance metric
    """
    # Codon count sums
    N_S = codon_count_sum(spike_dict)
    N_W = codon_count_sum(window_dict)

    # Compute distance
    distance = 0
    for codon in spike_dict.keys():
        F_S = spike_dict[codon] / N_S
        F_W = window_dict[codon] / N_W
        distance += (F_S - F_W)**2

    return distance


def codon_count_sum(codon_dict: dict) -> int:
    """Compute the total number codons observed in a sequence or window"""
    return sum(codon_dict.values())


def parse_FASTA_file(fasta_file: typing.IO) -> str:
    """Reads an input FASTA file to determine the genomic sequence.

    Considers that comment lines start with '>' and that sequences, if long 
    enough, can extend multiple lines.

    Parameters
    ----------
    fasta_fh : typing.IO
        Fasta file object

    Returns
    -------
    tuple[str, str]
        Pair of raw sequences to be processed
    """
    # Initialize data structure to hold sequence
    seq = []

    # Read FASTA line by line
    for line in fasta_file:
        if not line.startswith('>'):  # Comment line
            # append string to seq list
            seq.append(line.rstrip())

    # Merge strings into single string to create the sequence
    seq = ''.join(seq)
    return seq


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--spike_file', type=pathlib.Path,
                        help='Input file for the spike protein sequence')
    parser.add_argument('--viral_file', type=pathlib.Path,
                        help='Input file for the viral reference genome')
    parser.add_argument('--output_dir', type=pathlib.Path, default=pathlib.Path.cwd(),
                        help='The desired directory to create output files, by default the current working directory.')
    args = parser.parse_args()
    spike_file = args.spike_file
    viral_file = args.viral_file
    out_dir = args.output_dir

    # Sequences
    print('Loading sequences...')
    with open(spike_file) as f:
        spike_seq = parse_FASTA_file(f)
    with open(viral_file) as f:
        viral_seq = parse_FASTA_file(f)

    # Run functions
    print('Computing distances...')
    spike_dict = count_codon(spike_seq)
    window_distances = compute_window_distances(viral_seq, spike_dict)

    # Output spike codon count
    with open(out_dir.joinpath('spike_codon_count.txt'), 'w') as fh:
        for codon, count in sorted(spike_dict.items()):
            print(codon, count, sep='\t', file=fh)
    print(
        f'Exported spike codon counts to {out_dir.joinpath("spike_codon_count.txt")}')

    # Get positions and distances
    positions = []
    distances = []
    for pos, dist in window_distances.copy():
        positions.append(pos)
        distances.append(dist)

    # Print good_reads to file
    with open(out_dir.joinpath('window_distance_table.txt'), 'w') as fh:
        print('pos', 'val', sep='\t', file=fh)
        for pos_val in window_distances:
            print(*pos_val, sep='\t', file=fh)
    print(
        f'Exported distances to {out_dir.joinpath("window_distance_table.txt")}')

    # Plotting
    fig, ax = plt.subplots()
    ax.scatter(positions, distances, s=5, c='red', marker='.', linewidths=0)
    ax.set_xlabel('Window position')
    ax.set_ylabel('Difference distance')
    ax.set_title(
        'Difference distance between SARS-CoV-2 spike protein\nsequence and window at each starting position in full\nviral genome')
    plt.grid()
    plt.tight_layout()
    plt.savefig(out_dir.joinpath('window_difference.png'))
    print(f'Exported plot to {out_dir.joinpath("window_difference.png")}')
    # plt.show()
