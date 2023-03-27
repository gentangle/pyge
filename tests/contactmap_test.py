"""Test contactmap function."""
import csv
import pathlib

import numpy as np

import pyge.contacts.contactmap as cm


parent = pathlib.Path(__file__).parent.resolve()


with open(parent / "data/1srl_seqres_cut.csv") as in_file:
    reader = csv.reader(in_file)
    seq = list(reader)[0]

def compute_contactmap_sequence_test()
    """Test correctness contact map calculation"""
    expected = np.load(parent / "data/1ucs_cm.npy")
    out = cm.compute_contactmap(
        parent / "data/1ucs.pdb",
        1,
        "A",
        4.5,
        to_ignore=["HOH"]
    )
    assert np.allclose(out, expected, atol=0.0001)

def compute_contactmap_sequence_test():
    """Test computation of the contact map with sequence."""
    out = cm.compute_contactmap(
        parent / "data/1srl_w_missing_res.pdb",
        1,
        "A",
        4.5,
        to_ignore=["HOH"],
        sequence=seq,
    )

    assert out.shape == (64, 64)

    # test first rows
    assert not np.any(out[0])
    assert not np.any(out[1])
    assert not np.any(out[2])
    assert not np.any(out[3])
    assert not np.any(out[4])
    assert not np.any(out[5])
    assert not np.any(out[6])
    assert not np.any(out[7])
    # should be False, first res to be in the chain
    assert np.any(out[8])

    # residue removed in the middle (test on columns)
    assert np.any(out[:, 59])
    assert not np.any(out[:, 60])
    assert np.any(out[:, 59])
