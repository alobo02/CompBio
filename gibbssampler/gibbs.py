#!/usr/bin/env python3

"""Implements a simplified variant of the Gibbs sampling algorithm to find a 
common sequence motif.

TODO: finish module docstring

Command Line Usage
------------------
TODO: show command line usage
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
TC-CAAATAGAC
TCGCAAATATAC

$ python3 align.py -M 3 -m -2 -g -2 < tests/input/aligntest.input2 > tests/output/my.output2
$ python3 align.py --match_score 4 --gap_penalty -1 < tests/input/aligntest.input3 > tests/output/my.output3

Author: Alexander A. Lobo
"""

from ast import arg
import typing
import sys
import argparse
import numpy as np


class SimpleGibbsSampler:
    """An object-oriented instance of a Gibbs Sampler that is used to determine
    the MSA of the motif found from a set of sequences.

    Usage Examples
    --------------
    TODO: fix usage examples
    >>> align = GlobalSeqAlignment(x='ABC', y='AC', M=4, m=-2, g=-2)
    >>> align.get_optimal_alignment()
    ('ABC', 'A-C')
    """

    def __init__(
        self,
        seqs: list[str],
        motif_len: int = 6,
        lang: list[str] = ['A', 'C', 'G', 'T'],
        lang_prob: float = 0.25,
        seed: int = 20,
    ) -> None:
        """Instantiates an instance of the SimpleGibbsSampler class.

        Parameters
        ----------
        seqs : list[str]
            _description_
        motif_len : int, optional
            _description_, by default 6
        lang : list[str], optional
            _description_, by default ['A', 'C', 'G', 'T']
        lang_prob : float, optional
            _description_, by default 0.25
        seed : int, optional
            _description_, by default 20
        """
        # Arguments
        self.seqs = seqs
        self.motif_len = motif_len
        self.lang = lang
        self.lang_prob = lang_prob
        self.random_state = np.random.RandomState(seed)

        # Keep track of sequence positions positions
        self.positions = None

    def pick_init_positions(self) -> list[int]:
        """Generates and returns a list of random starting positions for each 
        sequence (i.e. set I = {i1, ..., ik}).

        Returns
        -------
        init_pos : list[int]
            List of random starting positions
        """

        init_pos = [
            self.random_state.randint(
                low=0,
                high=len(seq)-self.motif_len+1,
            ) for seq in self.seqs
        ]
        return init_pos

    def build_pseudo_count_matrix(self, motif) -> np.ndarray:
        """Given a list of string motifs, builds a counts matrix using 
        pseudocounts of 1 and returns it.

        Parameters
        ----------
        motif : _type_
            _description_

        Returns
        -------
        np.ndarray
            _description_
        """

        count_matrix = np.ones((len(self.lang), self.motif_len))
        for kmer in motif:
            for j in range(len(kmer)):
                i = self.lang.index(kmer[j])
                count_matrix[i, j] += 1
        return count_matrix

    def build_pssm(self, count_matrix) -> np.ndarray:
        """Given a counts matrix of type numpy array, builds a log-odds 
        position-specific scoring matrix (PSSM) and returns it.

        Parameters
        ----------
        count_matrix : _type_
            _description_

        Returns
        -------
        pssm : np.ndarray
            _description_
        """

        freq_matrix = count_matrix / \
            (count_matrix, np.sum(axis=0) + len(self.lang))
        pssm = np.log(freq_matrix / self.lang_prob)
        return pssm

    def score_seq_windows(self, seq, pssm) -> list[float]:
        """Calculates and returns a list of log-odds scores for each window 
        (i.e. possible motif) in the given sequence using given log-odds PSSM.

        Parameters
        ----------
        seq : _type_
            _description_
        pssm : _type_
            _description_

        Returns
        -------
        window_scores = list[float]
            _description_
        """
        window_scores = []
        for i in range(len(seq)-self.motif_len+1):
            kmer = seq[i:i+self.motif_len]
            window_score = self.get_window_score(kmer, pssm)
            window_scores.append(window_score)
        return window_scores

    def get_msa(self) -> list[str]:
        """Returns a list of string motifs according to the starting positions 
        configured for each sequence. 

        Must include all k many motifs for the k many sequences.

        Returns
        -------
        msa : list[str]
            _description_
        """
        msa = list(
            map(lambda seq, pos: seq[pos:pos+self.motif_len], self.seqs, self.positions))
        return msa

    def run_sampler(self) -> list[str]:
        """Runs all of the iterations of the sampler and returns a list of 
        string MSA.

        Returns
        -------
        list[str]
            _description_

        TODO
        ----
        * Need to figure out the optimal way to randomly choose s*
        * Don't think I can use choice because I also need to keep track of the positions
        * Think about how to avoid picking the same s*
        * To avoid accidently converging early, should keep track of another list
        that has "previous best positions"
        """
        same_position = False
        seqs_without_s_star = self.seqs.copy()
        while not same_position:
            # Pick I
            self.positions = self.pick_init_positions()

            # Randomly pick s*
            s_star = self.random_state.choice(seqs_without_s_star)

            # Remove s* from S
            seqs_without_s_star = self.seqs.copy().remove(s_star)

            # Build PSSM from S \ s*

            # Calculate the score for each window in s*

            # Choose starting position of s* as highest scoring

            # Compare best starting position in s* to previous time it was withheld

            # Iterate if not

    def get_window_score(self, kmer, pssm):
        score = 0
        for j, letter in enumerate(kmer):
            i = self.lang.index(letter)
            score += pssm[i, j]
        return score


if __name__ == "__main__":
    # TODO: fix the main code
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

    if not m > 2*g:
        raise argparse.ArgumentTypeError(
            'The mismatch score [-m, --mismatch_score] '
            'must be greater than 2 times the gap penalty [-g, --gap_penalty].',
        )

    # FASTA
    x, y = parse_FASTA_file(sys.stdin)

    # Class, method, and print
    align = GlobalSeqAlignment(x=x, y=y, M=M, m=m, g=g)
    print(*align.get_optimal_alignment(), sep='\n')
