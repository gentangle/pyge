"""Test the PDB parser module."""

import pathlib

import pytest

import pyge.contacts.pdb_parser as pdb_parser

parent = pathlib.Path(__file__).parent.resolve()


@pytest.mark.parametrize(
    "to_include, to_ignore, exp_res",
    [
        [None, ["HOH"], 64],
        [None, ["W"], 64],
        [None, ["HETATM"], 64],
        [["HOH"], None, 304],
        [["HOH"], ["HETATM"], 304],
    ],
)
def test_get_residues(to_include, to_ignore, exp_res):
    """Test the get_residues function."""
    res_list = pdb_parser.get_residues(
        str(parent / "data" / "1ucs.pdb"),
        1,
        "A",
        to_include=to_include,
        to_ignore=to_ignore,
        debug=False,
    )
    assert len(res_list) == exp_res
