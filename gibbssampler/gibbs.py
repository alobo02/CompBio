#!/usr/bin/env python3

"""Implements a simplified variant of the Gibbs sampling algorithm to find a 
common sequence motif.

Command Line Usage
------------------
The general format of running the module from the command line is as follows:

$ python3 gibbs.py [args] < [FASTA_input_file] > [output_file]

where

[args] can be any of the following arguments:
  -h, --help            show this help message and exit
  --motif_len MOTIF_LEN
                        The fixed motif length, default 6
  --lang LANG [LANG ...]
                        List characterizing the residue language symbols,
                        default ['A', 'C', 'G', 'T']
  --lang_prob LANG_PROB
                        Background probability of each residue, typically 1/len(lang),
                        default 0.25
  --seed SEED           Integer seed for random number generator, default 20

Run the following command for information on arguments:
$ python3 gibbs.py -h

[FASTA_input_file] is the file path of the input file used to import the
squences used to generate the MSA.

[output_file] is the desired file path for the print log

Usage Examples
--------------
$ python3 gibbs.py < tests/input/Gibbs.short.fasta > tests/output/my.output1
$ python3 gibbs.py --motif_len 5 --seed 14 < tests/input/Gibbs.fasta > tests/output/my.output2

Author: Alexander A. Lobo
"""

from ast import arg
import typing
import sys
import argparse
import numpy as np


class SimpleGibbsSampler:
    """An object-oriented instance of a Gibbs Sampler that is used to determine
    the MSA of the best motif found from a set of sequences.

    Usage Examples
    --------------
    >>> sgs = SimpleGibbsSampler(seqs=['AAA', 'AAC', 'ACG'], motif_len=2)
    >>> sgs.run_sampler()
    >>> sgs.get_msa()
    ['AA', 'AA', 'AC']

    TODO
    ----
    * fix usage examples
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
            List of sequences to be used to determine MSA of best motifs
        motif_len : int, optional
            Desired length of motif, by default 6
        lang : list[str], optional
            List characterizing the residue language symbols, by default ['A', 'C', 'G', 'T']
        lang_prob : float, optional
            Background probability of each residue occuring, by default 0.25
        seed : int, optional
            Integer value for random number generator, by default 20
        """
        # Arguments
        self.seqs = seqs
        self.motif_len = motif_len
        self.lang = lang
        self.lang_prob = lang_prob
        self.random_state = np.random.RandomState(seed)

        # Keep track of sequence positions positions, initialize with random gen
        self.positions = self.pick_init_positions()

    def pick_init_positions(self) -> list[int]:
        """Generates and returns a list of random starting positions for each 
        sequence (i.e. set I = {i1, ..., ik}).

        Returns
        -------
        init_pos : list[int]
            List of random starting positions
        """
        # Determine random starting positions of each sequence
        # Uses list comprehension to determine a random integer between 0 and
        # len(seq)-self.motif_len+1 for each seq
        init_pos = [
            self.random_state.randint(
                low=0,
                high=len(seq)-self.motif_len+1,
            ) for seq in self.seqs
        ]
        return init_pos

    def build_pseudo_count_matrix(self, motif: list[str]) -> np.ndarray:
        """Given a list of string motifs, builds a counts matrix using 
        pseudocounts of 1 and returns it.

        Parameters
        ----------
        motif : list[str]
            List of string motifs

        Returns
        -------
        count_matrix : np.ndarray
            Matrix representing counts of residues in each motif position, 
            dimensions: (len(lang) x motif_leng)
        """
        # Initialize count matrix with pseudo-counts=1
        count_matrix = np.ones((len(self.lang), self.motif_len))

        # Iterate through motif strings and add to count matrix accordingly
        # For each letter in kmer, determine index in lang and use to update matrix
        for kmer in motif:
            for j in range(len(kmer)):
                i = self.lang.index(kmer[j])
                count_matrix[i, j] += 1
        return count_matrix

    def build_pssm(self, count_matrix) -> np.ndarray:
        """Given a counts matrix, builds a log-odds position-specific scoring 
        matrix (PSSM) and returns it.

        Uses following formulas to fill frequency and PSSM matrices:
            f[u,b] = (n[u,b] + ε) / (n_tot + |Σ|*ε)
            pssm[u,b] = log2(f[u,b] / p[b])

        Parameters
        ----------
        count_matrix : np.ndarray
            Matrix representing counts of residues in each motif position, 
            dimensions: (len(lang) x motif_leng)

        Returns
        -------
        pssm : np.ndarray
            Position-specific scoring matrix that converts likelihood to log-odds
            scores using log-base2
        """
        freq_matrix = count_matrix / (len(self.seqs) - 1 + len(self.lang))
        pssm = np.log2(freq_matrix / self.lang_prob)
        return pssm

    def score_seq_windows(self, seq: str, pssm: np.ndarray) -> list[float]:
        """Calculates and returns a list of log-odds scores for each window 
        (i.e. possible motif) in the given sequence using given log-odds PSSM.

        Parameters
        ----------
        seq : str
            Sequence to score each window in using PSSM
        pssm : np.ndarray
            Position-specific scoring matrix

        Returns
        -------
        window_scores = list[float]
            List of scores for each window in seq
        """
        # Initialize empty list
        window_scores = []
        # Iterate through windows in seq and compute each score using helper
        for i in range(len(seq)-self.motif_len+1):
            kmer = seq[i:i+self.motif_len]
            window_score = self._get_window_score(kmer, pssm)
            window_scores.append(window_score)
        return window_scores

    def get_msa(self) -> list[str]:
        """Returns a list of string motifs according to the starting positions 
        configured for each sequence. 

        Must include all k many motifs for the k many sequences.

        Returns
        -------
        msa : list[str]
            List of strings that reflect current MSA
        """
        # Use the seq and positions lists to index the windows in each sequence
        # Uses anonymous function in map, indexing using elements in both lists
        msa = list(
            map(lambda seq, pos: seq[pos:pos+self.motif_len], self.seqs, self.positions))
        return msa

    def run_sampler(self) -> list[str]:
        """Runs all of the iterations of the sampler and returns a list of 
        string MSA.

        Terminating condition is when the best position determined in s* is the
        same position for that sequence that was determined previously.

        Returns
        -------
        best_msa : list[str]
            The MSA representing the best motif determined from a set of 
            sequences
        """
        # Intialize conditions for while loop
        same_position = False  # if s* position is the same as previous
        previous_s_star_ind = None  # The sequence number of the previous s*
        s_star_ind = None  # The sequence number of the current s*
        s_star_pos = None  # The index position in s* determined from PSSM scores

        # Set up an additional list that represents the previous positions
        # determined in s* for each sequence.
        # This avoid an early termination if s* pos happens to be the same as
        # the randomly chosen initial index.
        previous_positions = [None]*len(self.positions)

        # Set up iteration count
        iteration = 0

        # Iterate until position in s* is same as previous determined
        while not same_position:
            # Increment Iteration Number
            iteration += 1

            # Build MSA using every sequence
            msa = self.get_msa()

            # Randomly pick s* that is not the same as previous s*
            while previous_s_star_ind == s_star_ind:
                s_star_ind = self.random_state.choice(len(self.seqs))
            s_star = self.seqs.copy()[s_star_ind]

            # Remove s* from MSA
            msa_without_s_star = msa.copy()
            del msa_without_s_star[s_star_ind]

            # Build PSSM from S \ s*
            count_matrix = self.build_pseudo_count_matrix(msa_without_s_star)
            pssm = self.build_pssm(count_matrix)

            # Calculate the score for each window in s*
            window_scores = self.score_seq_windows(s_star, pssm)

            # Choose starting position of s* as highest scoring
            s_star_pos = np.argmax(window_scores)

            # Update previous s* index
            previous_s_star_ind = s_star_ind

            # Compare best starting position in s* to previous time it was withheld
            if not s_star_pos == previous_positions[s_star_ind]:
                # Update positions
                previous_positions[s_star_ind] = s_star_pos
                self.positions[s_star_ind] = s_star_pos
            else:
                same_position = True

            # Print log
            print(f'Iteration Number : {iteration}')
            print(f's* : {s_star}')
            print(f's* index from seqs (starts at 0) : {s_star_ind}')
            print('PSSM :')
            print(pssm)
            print('MSA :')
            print(*self.get_msa(), sep='\n', end='\n\n')

        # Obtain best MSA
        best_msa = self.get_msa()
        return best_msa

    def _get_window_score(self, kmer: str, pssm: np.ndarray) -> float:
        """Computes the score for a single window given a PSSM.

        Parameters
        ----------
        kmer : str
            String representing single window in sequence
        pssm : np.ndarray
            Position-specific scoring matrix

        Returns
        -------
        score : float
            Score of the window
        """
        score = 0
        for j, letter in enumerate(kmer):
            i = self.lang.index(letter)
            score += pssm[i, j]
        return score


def parse_FASTA_file(fasta_file: typing.IO) -> tuple[str, str]:
    """Reads an input FASTA file to determine the input sequences.

    Considers that comment lines start with '>' and that sequences, if long 
    enough, can extend multiple lines.

    Parameters
    ----------
    fasta_fh : typing.IO
        Fasta file object

    Returns
    -------
    list[str]
        List of raw sequences to be processed
    """
    # Initialize stack to hold sequences
    sequences = []

    # Read FASTA line by line
    for line in fasta_file:
        if line.startswith('>'):  # Comment line
            # Append to stack
            sequences.append([])
        else:
            # Append string to list at top of stack, removing newline characters
            sequences[-1].append(line.rstrip())

    # Merge strings into single string for each sequence and store in list
    sequences = list(map(lambda z: ''.join(z), sequences))
    return sequences


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--motif_len', type=int, default=6,
                        help='The fixed motif length')
    parser.add_argument('--lang', nargs='+', default=['A', 'C', 'G', 'T'],
                        help='List characterizing the residue language symbols')
    parser.add_argument('--lang_prob', type=float, default=0.25,
                        help='Background probability of each residue, typically 1/len(lang)')
    parser.add_argument('--seed', type=int, default=20,
                        help='Integer seed for random number generator')
    args = parser.parse_args()
    motif_len = args.motif_len
    lang = args.lang
    lang_prob = args.lang_prob
    seed = args.seed

    # FASTA
    seqs = parse_FASTA_file(sys.stdin)

    # Class and run method
    sgs = SimpleGibbsSampler(seqs, motif_len, lang, lang_prob, seed)
    msa = sgs.run_sampler()
