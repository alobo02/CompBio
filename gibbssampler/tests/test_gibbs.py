#!/usr/bin/env python3
import filecmp
from io import StringIO
import sys
import os
import pytest
import filecmp
import pathlib
from gibbssampler import gibbs


@pytest.mark.parametrize(
    "seed",
    [
        (20),
        (14),
        (1),
        (2),
        (3),
        (4),
        (5),
        (10),
        (15),
        (30),
    ],
)
def test_module(seed):
    module_path = pathlib.Path(__file__).parents[1].joinpath('gibbs.py')
    input_path = pathlib.Path(__file__).parent.joinpath(
        'input', 'Gibbs.short.fasta')
    test_output_path = pathlib.Path(__file__).parent.joinpath(
        'output', f'expected_output.seed{seed}')
    my_output_path = pathlib.Path(__file__).parent.joinpath(
        'output', f'output.seed{seed}')
    os.system(
        f"python3 {module_path} --seed {seed} < {str(input_path)} > {str(my_output_path)}")
    assert filecmp.cmp(test_output_path, my_output_path)
