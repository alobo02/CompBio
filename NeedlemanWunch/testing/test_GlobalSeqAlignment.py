#!/usr/bin/env python3
import filecmp
from io import StringIO
import sys
import os
import pytest
import filecmp
import pathlib


@pytest.mark.parametrize(
    "test_case",
    [
        ("1"),
        ("2"),
        ("3"),
    ],
)
def test_align(test_case):
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
