#!/usr/bin/env python3

"""Implements a global DNA sequence alignment using a simple similarity scoring 
scheme with a linear gap penalty: 
* M for each matched pair in the alignment
* m for each mismatch
* g for each gap

The module can either be imported as a package to use the
GlobalSequenceAlignment class and class methods or run directly from the command
line.

Command Line Usage
------------------
The general format of running the module from the command line is as follows:

$ python3 align.py [args] < [FASTA_input_file] > [output_file]

where

[args] can be any of the following arguments:
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

[FASTA_input_file] is the file path of the input file used to import the
squences to be aligned

[output_file] is the desired file path of the file to print the optimal global
aligments to

Usage Examples
--------------
$ python3 align.py < tests/input/aligntest.input1 > tests/output/my.output1
$ python3 align.py -M 3 -m -2 -g -2 < tests/input/aligntest.input2 > tests/output/my.output2
$ python3 align.py --match_score 4 --gap_penalty -1 < tests/input/aligntest.input3 > tests/output/my.output3

Author: Alexander A. Lobo
"""

from ast import arg
import typing
import sys
import argparse
import numpy as np


class GlobalSequenceAlignment:
    """An object-oriented instance of an alignment than can be used to find
    the optimal global alignment of two sequences.

    Usage Examples
    --------------
    >>> align = GlobalSequenceAlignment(x='ABC', y='AC', M=4, m=-2, g=-2)
    >>> align.get_optimal_alignment()
    ('ABC', 'A-C')
    """

    def __init__(
        self,
        x: str,
        y: str,
        M: float = 4,
        m: float = -2,
        g: float = -2
    ) -> None:
        """Instantiates an instance of the GlobalSequenceAlignment class.

        Parameters
        ----------
        x : str
            First sequence
        y : str
            Second sequence
        M : float, optional
            Match score, by default 4
        m : float, optional
            Mismatch penalty, by default -2
        g : float, optional
            Linear gap penalty, by default -2
        """

        # Strings
        self.x = '_' + x
        self.y = '_' + y

        # Scoring Scheme
        self.M = M
        self.m = m
        self.g = g

        # Obtain dimension of matrix
        m_dim, n_dim = len(self.x), len(self.y)
        self.m_dim = m_dim
        self.n_dim = n_dim

        # Initialize alignment (A) matrix
        A = np.zeros((m_dim, n_dim))
        A[:, 0] = np.arange(0, g*m_dim, g) # First column gap scores
        A[0, :] = np.arange(0, g*n_dim, g) # First row gap scores
        self.A = A

        # Initialize pointer (P) matrix
        P = np.zeros((m_dim, n_dim))
        P[1:, 0] = 2 # First column always points "up"
        P[0, 1:] = 3 # First row always points "left"
        self.P = P

        # Compute matrix scores
        self._compute_matrix()

    def _compute_matrix(self) -> None:
        """Fills in the alignment matrix (A).

        Iterates over each index in A and updated the score.
        """
        for i in range(1, self.m_dim):
            for j in range(1, self.n_dim):
                self._update_matrix(i, j)

    def _update_matrix(self, i: int, j: int) -> None:
        """Updates the score of the alignment matrix at A(i,j) using the
        scoring scheme.

        Also updated the pointer matrix at P(i,j) to indicate which direction
        was used to compute the score

        Parameters
        ----------
        i : int
            Alignment matrix row index
        j : int
            Alignment matrix column index
        """
        cases = [
            self.A[i-1, j-1] + self._score(self.x[i], self.y[j]),
            self.A[i-1, j] + self.g,
            self.A[i, j-1] + self.g
        ]
        self.A[i, j] = max(cases)
        self.P[i, j] = np.argmax(cases) + 1

    def _score(self, x_i: str, y_j: str) -> float:
        """Determines the score of two characters in the alphabet (excluding
        gaps)

        Parameters
        ----------
        x_i : str
            First character
        y_j : str
            Second character

        Returns
        -------
        float
            Score, M or m, depending on whether the characters match or not
        """
        if x_i == y_j:
            return self.M
        else:
            return self.m

    def get_optimal_alignment(self) -> tuple[str, str]:
        """Traces back through the pointer matrix (P) to determine the optimal
        global alignment.

        Builds a stack of both sequences in reverse order until either string
        index reaches zero.

        Returns
        -------
        tuple[str, str]
            Optimal alignment pair of both sequences
        """
        i = self.m_dim-1
        j = self.n_dim-1

        x_stack = []
        y_stack = []

        while i > 0 or j > 0:
            pointer = self.P[i, j]
            if pointer == 1:
                x_stack.append(self.x[i])
                y_stack.append(self.y[j])
                i -= 1
                j -= 1
            elif pointer == 2:
                x_stack.append(self.x[i])
                y_stack.append('-')
                i -= 1
            elif pointer == 3:
                x_stack.append('-')
                y_stack.append(self.y[j])
                j -= 1

        x_aligned = ''.join(x_stack[::-1])
        y_aligned = ''.join(y_stack[::-1])

        return (x_aligned, y_aligned)


def parse_FASTA_file(fasta_fh: typing.IO) -> tuple[str, str]:
    """Reads an input FASTA file to determine the input sequences to be aligned.

    Parameters
    ----------
    fasta_fh : typing.IO
        Fasta file object

    Returns
    -------
    tuple[str, str]
        Pair of raw sequences to be processed
    """
    x_y = [[], []]
    i = -1
    for line in fasta_fh.readlines():
        if line.startswith('>'):
            i += 1
        else:
            x_y[i].append(line.rstrip())
    x_y = tuple(map(lambda z: ''.join(z), x_y))
    return x_y


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-M', '--match_score', type=float, default=4,
                        help='Match score to be used for global alignment')
    parser.add_argument('-m', '--mismatch_score', type=float, default=-2,
                        help='Mismatch score to be used for global alignment')
    parser.add_argument('-g', '--gap_penalty', type=float, default=-2,
                        help='Linear gap penalty to be used for global alignment')
    args = parser.parse_args()
    M = args.match_score
    m = args.mismatch_score
    g = args.gap_penalty

    # FASTA
    x, y = parse_FASTA_file(sys.stdin)

    # Class & Method
    align = GlobalSequenceAlignment(x=x, y=y, M=M, m=m, g=g)
    print(*align.get_optimal_alignment(), sep='\n')
