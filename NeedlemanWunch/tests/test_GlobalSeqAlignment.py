#!/usr/bin/env python3
import filecmp
from io import StringIO
import sys
import os
import pytest
import filecmp
import pathlib
from NeedlemanWunch import align


@pytest.mark.parametrize(
    "x, y, expected",
    [
        ('ABC', 'ABC', ('ABC', 'ABC')),
        ('ABC', 'AB', ('ABC', 'AB-')),
        ('ABC', 'A', ('ABC', 'A--')),
        ('ABC', '', ('ABC', '---')),
        ('', '', ('', '')),
        ('A', 'ABC', ('A--', 'ABC')),
        ('ABCDEF', 'AF', ('ABCDEF', 'A----F')),
        ('ABCDEF', 'A', ('ABCDEF', 'A-----')),
        ('ABCDEF', 'F', ('ABCDEF', '-----F')),
        ('ABCDEF', 'ACF', ('ABCDEF', 'A-C--F')),
        ('ABCDEF', 'QCWERTY', ('ABCDE--F', '-QCWERTY')),
    ],
)
def test_get_optimal_alignment(x, y, expected):
    a_align = align.GlobalSequenceAlignment(x=x, y=y)
    assert a_align.get_optimal_alignment() == expected


@pytest.mark.parametrize(
    "test_case",
    [
        ("1"),
        ("2"),
        ("3"),
    ],
)
def test_module(test_case):
    module_path = pathlib.Path(__file__).parents[1].joinpath('align.py')
    input_path = pathlib.Path(__file__).parent.joinpath(
        'input', f'aligntest.input{test_case}')
    test_output_path = pathlib.Path(__file__).parent.joinpath(
        'output', f'aligntest.output{test_case}')
    my_output_path = pathlib.Path(__file__).parent.joinpath(
        'output', f'my.output{test_case}')
    os.system(
        f"python3 {module_path} < {str(input_path)} > {str(my_output_path)}")
    assert filecmp.cmp(test_output_path, my_output_path)
