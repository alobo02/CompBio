#!/usr/bin/env python3
"""Partial assembly of a set of DNA sequence reads.

Implements a de Bruijn graph for error detection and assembly to be used on 
a provided set of sequence reads.

Command Line Usage
------------------
The general format of running the module from the command line is as follows:

$ python3 assembler.py [args] < [reads_input_file]

where

[args] can be any of the following arguments:
  -h, --help            show help message and exit
  --k K                 The length of each k-mer
  --min_length MIN_LENGTH
                        Minimum length desired of contigs to be returned
  --output_dir OUTPUT_DIR
                        The desired directory to create output files, by default the current
                        working directory.

Run the following command for information on arguments:
$ python3 assembler.py -h

[reads_input_file] is the file path of the input file used to import the
squences reads to assemble the DNA sequence.

Usage Examples
--------------
$ python3 assembler.py < tests/input/sequence_reads
$ python3 assembler.py --k 31 --min_length 100 < tests/input/sequence_reads

Author: Alexander A. Lobo
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
    """An object-oriented instance of a de Bruijn Graph.

    Uses a dictionary to keep track of the nodes of the graph, where each
    k-mer is a key, and an instance of a GraphNode object is the corresponding 
    value.

    Usage Examples
    --------------
    >>> graph = DeBruijnGraph(['AAABCDAAABCD', 'BCDAB', 'ABC'], k=3)
    >>> graph.build_graph()
    >>> graph.get_contigs()
    (['DAAABCD', 'DAB'], [7, 3])
    """

    def __init__(
        self,
        reads: list[str],
        k: int = 31,
    ) -> None:
        """Instantiates an instance of the de Bruijn Graph class.

        Parameters
        ----------
        reads : list[str]
            List of sequence reads
        k : int, optional
            Desired K-mer size, by default 31
        """
        # Parameters
        self.reads = reads
        self.k = k

        # Set up an empty dictionary to store the the nodes
        self.nodes = dict()

    def build_graph(self) -> None:
        """Builds the de Bruijn graph out of the given list of reads by looping
        through each read and looping through each read's k-mers to add nodes 
        and edges.
        """
        # Loop through reads
        for read in self.reads:
            # Loop through kmers
            kmers = self._get_kmers(read)
            for i, kmer in enumerate(kmers):
                # Add each kmer as a node
                self.add_node(kmer)
                # If on the 2nd kmer or after, add the in/out edges to the
                # previous kmer and current kmer node pair
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
        """
        # Initialize empty list of kmers
        kmers = []
        # Detmine the number of kmers in the read
        num_kmers = len(read) - self.k + 1
        # Iterate through kmers and append to list
        for i in range(num_kmers):
            kmer = read[i:i+self.k]
            kmers.append(kmer)
        return kmers

    def add_node(self, kmer: str) -> None:
        """Adds the string k-mer value, kmer, to the graph.

        Creates a new node if it doesn't already exists, otherwise will update
        the count of the node by adding one to it.

        Parameters
        ----------
        kmer : str
            kmer to be added to the graph
        """
        if kmer not in self.nodes:
            self.nodes[kmer] = GraphNode(kmer)
        else:
            self.nodes[kmer].count += 1

    def add_in_edge(self, kmer: str, in_kmer: str) -> None:
        """Adds in_kmer as kmer node's incoming edge by appending in_kmer to its
        incoming list.

        Updates branch status if node is a branch after the update.

        Parameters
        ----------
        kmer : str
            Key of the desired node incoming edge list to update
        in_kmer : str
            Key of node to use to update kmer incoming edge list with
        """
        # Check to make sure that the in_kmer is not already an incoming edge
        if in_kmer not in self.nodes[kmer].incoming:
            # Update kmer incoming edge list
            self.nodes[kmer].incoming.append(in_kmer)
        # Update the node's branch status
        self._update_branch_node_status(kmer)

    def add_out_edge(self, kmer: str, out_kmer: str) -> None:
        """Adds out_kmer as kmer node's outgoing edge by appending out_kmer to
        its outgoing list.

        Updates branch status if node is a branch after the update.

        Parameters
        ----------
        kmer : str
            Key of the desired node outgoing edge list to update
        out_kmer : str
            Key of node to use to update kmer outgoing edge list with
        """
        # Check to make sure that the in_kmer is not already an outgoing edge
        if out_kmer not in self.nodes[kmer].outgoing:
            # Update kmer outgoing edge list
            self.nodes[kmer].outgoing.append(out_kmer)
        # Update the node's branch status
        self._update_branch_node_status(kmer)

    def remove_node(self, kmer: str) -> None:
        """Removes the node that represents kmer from the graph along with any
        edges to and from the removed node.

        Parameters
        ----------
        kmer : str
            Key of the desired node to remove
        """
        # Delete the dictionary item (assumes node exists)
        del self.nodes[kmer]
        # Iterate through rest of nodes and remove any incoming/outgoing edges
        # that have kmer
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
        """
        # Initialize empty list to store good reads
        good_reads = []
        # Iterate through reads
        for read in self.reads:
            # Get the kmers of the read
            kmers = self._get_kmers(read)
            # Assume that all kmers in the read are "good" (occurs more than once)
            kmers_all_good = True
            # Keep track of kmer index
            i = 0
            # Iterate through all kmers in read unless a 1 count is found
            while kmers_all_good and i < len(kmers):
                # Stop iterating if kmer with count = 1 is found
                if self.nodes[kmers[i]].count == 1:
                    kmers_all_good = False
                else:
                    i += 1
            # Include the read as a "good read" if all the kmers had count > 1
            if kmers_all_good:
                good_reads.append(read)
        return good_reads

    def get_contigs(self, min_length=0) -> None:
        """Finds a list of contigs that is greater than or equal to min_len 
        using the graph.

        Parameters
        ----------
        min_length : int, optional
            Desired minimum lengths of contigs, by default 0
        """
        # Remove branch nodes
        self._remove_branch_nodes()
        # Obtain the source nodes
        source_nodes = self._get_source_nodes()
        # For each source node, traverse the graph' and obtain the contigs
        base_contigs = [self._traverse_contig_from(source_node)
                        for source_node in source_nodes]
        # Filter contigs that are greater than or equal to min_length
        contigs = [contig for contig in base_contigs if len(
            contig) >= min_length]
        # Obtain the lengths of each of the contigs
        contig_lengths = list(map(lambda z: len(z), contigs))
        return contigs, contig_lengths

    def _remove_branch_nodes(self) -> None:
        """Removes all the branch nodes from the graph.

        The graph can now be thought of as graph'.
        """
        # Iterate through nodes in graph
        # Makes a copy of the keys to iterate through since dictionary will change
        for kmer in list(self.nodes.keys()):
            # If the node is a branch node, remove it
            if self.nodes[kmer].is_branch:
                self.remove_node(kmer)

    def _get_source_nodes(self) -> list[str]:
        """Obtains the source nodes in the graph.

        Source node have no incoming edges.

        Returns
        -------
        source_nodes : list[str]
            List of keys of source nodes
        """
        # Iterate through node dictionary and return keys that have empty
        # incoming edge lists
        source_nodes = [k for k, v in self.nodes.items() if not v.incoming]
        return source_nodes

    def _traverse_contig_from(self, source_kmer: str) -> str:
        """Traverse through the graph from a source node to onbtain the contig.

        Will iterate through graph by looking at outgoing edges until a sink is
        found. Sink nodes have no outgoing edges.

        Parameters
        ----------
        source_kmer : str
            Key of the source node

        Returns
        -------
        contig : str
            Contig obtained from traversing the graph
        """
        # Initialize contig with the source node kmer
        contig = source_kmer
        # Obtain the source node
        node = self.nodes[source_kmer]
        # Iterate through graph until a sink node is found
        while node.outgoing:
            # Next kmer is the outgoing edge of current node
            kmer = node.outgoing[-1]
            # Build the contig with last residue of next kmer
            contig += kmer[-1]
            # Update the node with the next contig
            node = self.nodes[kmer]
        return contig

    def _update_branch_node_status(self, kmer: str) -> None:
        """Checks the incoming and outgoing edges and updates the branch status
        of the node if either have more than 1 edge.

        Parameters
        ----------
        kmer : str
            Key of desired node to check

        Raises
        ------
        KeyError
            Raises an exception if trying to update the states of a node that 
            doesn't exist
        """
        if kmer not in self.nodes:
            raise KeyError(
                f'The provided k-mer "{kmer}" does not exist as a node in the graph.')
        # Update branch status if incoming/outgoing length > 1
        if len(self.nodes[kmer].incoming) > 1 or len(self.nodes[kmer].outgoing) > 1:
            self.nodes[kmer].is_branch = True


def parse_reads_file(reads_file: Union[typing.IO, list]) -> list[str]:
    """Reads an input file to determine the sequence reads.

    Assumes each new line is a sequence read.

    Parameters
    ----------
    reads_file : Union[typing.IO, list]
        File handle of the file containing sequence reads
        or a list containing the sequence reads

    Returns
    -------
    reads : list[str]
        List of the sequence reads
    """
    # Initialize stack to hold sequence reads
    reads = []

    # Read file line by line
    for line in reads_file:
        reads.append(line.rstrip())

    return reads


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=31,
                        help='The length of each k-mer')
    parser.add_argument('--min_length', type=int, default=0,
                        help='Minimum length desired of contigs to be returned')
    parser.add_argument('--output_dir', type=pathlib.Path, default=pathlib.Path.cwd(),
                        help='The desired directory to create output files, by default the current working directory.')
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
    print('Filtering good reads...')
    good_reads = graph.get_good_reads()

    # Print good_reads to file
    with open(out_dir.joinpath('good_reads'), 'w') as fh:
        print(*good_reads, sep='\n', file=fh)
    print(f'Exported good reads to {out_dir.joinpath("good_reads")}')

    # Good reads graph
    print('Building graph from good reads...')
    good_graph = DeBruijnGraph(good_reads, k)
    good_graph.build_graph()
    contigs, contig_lengths = good_graph.get_contigs(min_length=min_length)

    # Export contigs and contig lengths
    with open(out_dir.joinpath('output_contigs'), 'w') as fh:
        # Sort the contigs in the same order as the contig lengths
        print(*[contig for _, contig in sorted(zip(contig_lengths,
              contigs))], sep='\n', file=fh)
    with open(out_dir.joinpath('contig_lengths'), 'w') as fh:
        # Sort the contig lengths
        print(*sorted(contig_lengths), sep='\n', file=fh)
    print(f'Exported contigs to {out_dir.joinpath("output_contigs")}')
    print(f'Exported contig lengths to {out_dir.joinpath("contig_lengths")}')
