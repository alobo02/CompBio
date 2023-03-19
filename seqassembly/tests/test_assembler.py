#!/usr/bin/env python3
import filecmp
from io import StringIO
import sys
import os
import pytest
import filecmp
import pathlib
from seqassembly import assembler


@pytest.mark.parametrize(
    "testcase, k, min_length",
    [
        ('demo', 5, 8),
        ('main', 31, 100),
    ],
)
def test_module(testcase, k, min_length):
    module_path = pathlib.Path(__file__).parents[1].joinpath('assembler.py')
    input_path = pathlib.Path(__file__).parent.joinpath(
        'input', testcase, 'sequence_reads')
    output_path = pathlib.Path(__file__).parent.joinpath(
        'output', testcase)
    os.system(
        f"python3 {module_path} --k {k} --min_length {min_length} --out_dir {str(output_path)} < {str(input_path)}")
    assert filecmp.cmp(output_path.joinpath(
        'good_reads_expected'), output_path.joinpath('good_reads'))
    assert filecmp.cmp(output_path.joinpath(
        'output_contigs_expected'), output_path.joinpath('output_contigs'))
    assert filecmp.cmp(output_path.joinpath(
        'contig_lengths_expected'), output_path.joinpath('contig_lengths'))


@pytest.mark.parametrize(
    "k, min_length",
    [
        (31, 0),
    ],
)
def test_written_questions(k, min_length):
    module_path = pathlib.Path(__file__).parents[1].joinpath('assembler.py')
    input_path = pathlib.Path(__file__).parent.joinpath(
        'input', 'main', 'sequence_reads')
    output_path = pathlib.Path(__file__).parent.joinpath(
        'output', 'written_questions')
    os.system(
        f"python3 {module_path} --k {k} --min_length {min_length} --out_dir {str(output_path)} < {str(input_path)}")
