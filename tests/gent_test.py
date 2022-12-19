"""Test for the gent module"""
import pathlib
import json

import numpy as np
import MDAnalysis as mda

from pyge import gent

parent = pathlib.Path(__file__).parent.resolve()

ge_1ucs_native = np.load(
    parent / "data/GE_1ucs_native.npy",
    allow_pickle=True)
with open(parent / "data/1ucs_ge_configurations_modes.json", "r", encoding="utf-8") as fin:
    ge_1ucs_native_modes = json.load(fin)
cm_1ucs = np.load(
    parent / "data/1ucs_cm.npy"
)
u_pdb = mda.Universe(str(parent / "data/1ucs.pdb"))
ca_positions = u_pdb.select_atoms('name CA').positions


def test_gaussian_entanglement():
    """
    Testing the gaussian entanglement calculation for all possible loops
    """
    ge_computed = gent.ge_loops(cm_1ucs, ca_positions, 10, backend='cython')

    for ge_c, ge_true in zip(ge_computed, ge_1ucs_native):
        assert np.isclose(
            ge_c[2],
            ge_true[2]
        ), print(ge_c, ge_true)

def test_ge_configuration():
    """
    Test the possible selections for the GE of the whole configuration
    """
    ge_loops_computed = gent.ge_loops(cm_1ucs, ca_positions, 10, backend='cython')
    modes = ("max", "average")

    for mode in modes:
        ge_selected = gent.ge_configuration(ge_loops_computed, 10, mode)
        assert np.isclose(ge_selected[2], ge_1ucs_native_modes[mode], atol=0.0001)
        if mode == 'max':
            assert ge_selected[0][0] == ge_1ucs_native_modes["max_loop"][0]
            assert ge_selected[0][1] == ge_1ucs_native_modes["max_loop"][1]
            assert ge_selected[1][0] == ge_1ucs_native_modes["max_thr"][0]
            assert ge_selected[1][1] == ge_1ucs_native_modes["max_thr"][1]
