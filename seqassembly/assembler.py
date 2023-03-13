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

TODO
----
* Update module docstrring
"""

from ast import arg
import typing
from typing import Union
import sys
import argparse
import numpy as np
import pathlib


class GraphNode:
    """An object-oriented instance of a graph node in a de Bruijn Graph.

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
            value: str,
            incoming: list[str] = None,
            outgoing: list[str] = None,
            is_branch: bool = False,
            is_visited: bool = False,
            count: int = 1,
    ) -> None:
        """Instantiates an instance of the GraphNode class.

        Parameters
        ----------
        value : str
            K-mer value as string
        incoming : list[str], optional
            Incoming list of nodes connected to this node, by default None
        outcoming : list[str], optional
            Outgoing list of nodes connected from this node, by default None
        is_branch : bool, optional
            Whether the node is a branch node, by default False
        is_visited : bool, optional
            Whether the node has been visited when traversing for contigs, 
            by default False
        count : int, optional
            Number of times the node has been encountered, by default 1

        TODO
        ----
        * update docstring
        """
        # Arguments
        self.value = value
        if incoming is None:
            self.incoming = []
        else:
            self.incoming = incoming
        if outgoing is None:
            self.outgoing = []
        else:
            self.outgoing = outgoing
        self.is_branch = is_branch
        self.is_visited = is_visited
        self.count = count


class DeBruijnGraph:
    """_summary_
    """

    def __init__(
        self,
        reads: list[str],
        k: int = 31,
    ) -> None:
        # Parameters
        self.reads = reads
        self.k = k
        self.nodes = dict()

    def build_graph(self) -> None:
        """Builds the de Bruijn graph out of the given list of reads by looping
        through each read and looping through each read's k-mers to add nodes 
        and edges.

        TODO
        ----
        * Loop through reads
        * Loop through kmers in reads
        * Add each kmer as a node
        * Insert out edge for each kmer using next kmer
        * Insert in edge for next kmer using current kmer
        """
        # Loop through reads
        for read in self.reads:
            # Loop through kmers
            kmers = self._get_kmers(read)
            for i, kmer in enumerate(kmers):
                self.add_node(kmer)
                if i >= 1:
                    self.add_in_edge(kmer, kmers[i-1])
                    self.add_out_edge(kmers[i-1], kmer)

    def _get_kmers(self, read: str) -> list[str]:
        """Obtains a list of k-mers from a given read.

        Parameters
        ----------
        read : str
            A string representation of a read.

        Returns
        -------
        kmers : list[str]
            List of k-mers

        TODO
        ----
        * test
        * comment
        """
        kmers = []
        num_kmers = len(read) - self.k + 1
        for i in range(num_kmers):
            kmer = read[i:i+self.k]
            kmers.append(kmer)
        return kmers

    def add_node(self, kmer: str) -> None:
        """Adds the string k-mer value, kmer, to the graph.

        Creates a new node if it doesn't already exists, otherwise will update
        the count of the node.

        Parameters
        ----------
        kmer : str
            kmer to be added to the graph

        TODO
        ----
        * Implement
        """
        if kmer not in self.nodes:
            self.nodes[kmer] = GraphNode(kmer)
        else:
            self.nodes[kmer].count += 1

    def add_in_edge(self, kmer: str, in_kmer: str) -> None:
        """Adds in_kmer as kmer node's incoming edge by appending in_kmer to its
        incoming list.

        Parameters
        ----------
        kmer : str
            _description_
        in_kmer : str
            _description_

        TODO
        ----
        * Implement
        * Change is_branch to True if len(incoming) > 1
        """
        if in_kmer not in self.nodes[kmer].incoming:
            self.nodes[kmer].incoming.append(in_kmer)
        self._update_branch_node_status(kmer)

    def add_out_edge(self, kmer: str, out_kmer: str) -> None:
        """Adds out_kmer as kmer node's outgoing edge by appending out_kmer to
        its outgoing list.

        Parameters
        ----------
        kmer : str
            _description_
        out_kmer : str
            _description_

        TODO
        ----
        * Implement
        * Change is_branch to True if len(outgoing) > 1
        """
        if out_kmer not in self.nodes[kmer].outgoing:
            self.nodes[kmer].outgoing.append(out_kmer)
        self._update_branch_node_status(kmer)

    def remove_node(self, kmer: str) -> None:
        """Removes the node that represents kmer from the graph along with any
        edges to and from the removed node.

        Parameters
        ----------
        kmer : str
            _description_

        TODO
        ----
        * Remove node from dictionary
        * Iterate through all nodes and remove in/out edges that include kmer
          using list.remove()
        """
        del self.nodes[kmer]
        for k in self.nodes:
            if kmer in self.nodes[k].incoming:
                self.nodes[k].incoming.remove(kmer)
            if kmer in self.nodes[k].outgoing:
                self.nodes[k].outgoing.remove(kmer)

    def get_good_reads(self) -> list[str]:
        """Uses the graph to filter out reads that obtain kmers that occur only
        once and returns all the good reads that are left.

        Returns
        -------
        good_reads : list[str]
            A list of the good reads.

        TODO
        ----
        * iterate through reads
        * use while loop to iterate throught the kmers in the read
        * if a kmer with count = 1 is found, mark the read as bad and do not
          append to good_read list
        """
        good_reads = []
        for read in self.reads:
            kmers = self._get_kmers(read)
            kmers_all_good = True
            i = 0
            while kmers_all_good and i < len(kmers):
                if self.nodes[kmers[i]].count == 1:
                    kmers_all_good = False
                else:
                    i += 1
            if kmers_all_good:
                good_reads.append(read)
        return good_reads

    def get_contigs(self, min_length=0) -> None:
        """Finds a list of contigs that is equal to in length or longer than 
        min_len using the graph.

        Parameters
        ----------
        min_length : int, optional
            _description_, by default 0

        TODO
        ----
        * helper functions:
            * get_branch_nodes() or remove_branch_nodes()
            * get_source_nodes()
            * traverse_contig_from(source_kmer)
        * Need to remove branch nodes
        * Find sources
        * Traverse graph from source to sinks
        * Build contigs during traversal
        """
        self._remove_branch_nodes()
        source_nodes = self._get_source_nodes()
        base_contigs = [self._traverse_contig_from(
            source_node) for source_node in source_nodes]
        # what is this in length that is mentioned?
        contigs = [contig for contig in base_contigs if len(
            contig) >= min_length]
        contig_lengths = list(map(lambda z: len(z), contigs))
        return contigs, contig_lengths

    def _remove_branch_nodes(self) -> None:
        for kmer in list(self.nodes.keys()):
            if self.nodes[kmer].is_branch:
                self.remove_node(kmer)

    def _get_source_nodes(self) -> list[str]:
        source_nodes = [k for k, v in self.nodes.items() if not v.incoming]
        return source_nodes

    def _traverse_contig_from(self, source_kmer) -> str:
        contig = source_kmer
        node = self.nodes[source_kmer]
        while node.outgoing:
            kmer = node.outgoing[-1]
            contig += kmer[-1]
            node = self.nodes[kmer]
        return contig

    def _update_branch_node_status(self, kmer) -> None:
        if kmer not in self.nodes:
            raise KeyError(
                f'The provided k-mer "{kmer}" does not exist as a node in the graph.')
        if len(self.nodes[kmer].incoming) > 1 or len(self.nodes[kmer].outgoing) > 1:
            self.nodes[kmer].is_branch = True


def parse_reads_file(reads_file: Union[typing.IO, list]) -> list[str]:
    """Reads an input file to determine the sequence reads.

    Considers that comment lines start with '>' and that sequences, if long 
    enough, can extend multiple lines.

    Parameters
    ----------
    reads_file : Union[typing.IO, list]
        File handle of the file containing sequence reads
        or a list containing the sequence reads

    Returns
    -------
    reads : list[str]
        List of the sequence reads

    TODO
    ----
    * Test to make sure that it's working with the input format
    """
    # Initialize stack to hold sequence reads
    reads = []

    # Read FASTA line by line
    for line in reads_file:
        # if line.startswith('>'):  # Comment line
        #     # Append to stack
        #     reads.append([])
        # else:
        #     # Append string to list at top of stack, removing newline characters
        #     reads[-1].append(line.rstrip())
        reads.append(line.rstrip())

    # Merge strings into single string for each sequence and store in list
    # reads = list(map(lambda z: ''.join(z), reads))
    return reads


if __name__ == "__main__":
    # TODO:
    # * update run if module run directlt

    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=31,
                        help='The length of each k-mer')
    parser.add_argument('--min_length', type=int, default=0,
                        help='Minimum length desired of contigs to be returned')
    parser.add_argument('--output_dir', type=pathlib.Path, default=pathlib.Path.cwd(),
                        help='The desired directory to create outpur files, by default the current working directory.')
    args = parser.parse_args()
    k = args.k
    min_length = args.min_length
    out_dir = args.output_dir

    # Reads
    print('Obtaining reads...')
    reads = parse_reads_file(sys.stdin)

    # Instantiate graph from reads
    print('Building graph...')
    graph = DeBruijnGraph(reads, k)
    graph.build_graph()
    print('Obtaining good reads...')
    good_reads = graph.get_good_reads()

    # Print good_reads to file
    with open(out_dir.joinpath('good_reads'), 'w') as fh:
        print(*good_reads, sep='\n', file=fh)
    print(f'Exported good reads to {out_dir.joinpath("good_reads")}')

    # Good reads graph
    print('Building good graph...')
    good_graph = DeBruijnGraph(good_reads, k)
    good_graph.build_graph()
    contigs, contig_lengths = good_graph.get_contigs(min_length=min_length)

    # Export contigs and contig lengths
    with open(out_dir.joinpath('output_contigs'), 'w') as fh:
        print(*contigs, sep='\n', file=fh)
    with open(out_dir.joinpath('contig_lengths'), 'w') as fh:
        print(*sorted(contig_lengths), sep='\n', file=fh)
    print('Exported contigs and contig lengths')
