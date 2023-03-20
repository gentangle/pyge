"""Test for the gent module."""
import json
import logging
import pathlib
import sys

import MDAnalysis as mda
import numpy as np

from pyge import gent


logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

parent = pathlib.Path(__file__).parent.resolve()

ge_1ucs_native = np.load(parent / "data/GE_1ucs_native.npy", allow_pickle=True)
with open(
    parent / "data/1ucs_ge_configurations_modes.json", "r", encoding="utf-8"
) as fin:
    ge_1ucs_native_modes = json.load(fin)
cm_1ucs = np.load(parent / "data/1ucs_cm.npy")
u_pdb = mda.Universe(str(parent / "data/1ucs.pdb"))
ca_positions = u_pdb.select_atoms("name CA").positions


def test_gaussian_entanglement():
    """Testing the gaussian entanglement calculation for all possible loops."""
    ge_computed = gent.ge_loops(cm_1ucs, ca_positions, 10, backend="cython")

    for ge_comp, ge_exp in zip(ge_computed, ge_1ucs_native):
        if abs(ge_comp.n_term.value) >= abs(ge_comp.c_term.value):
            ge2com = ge_comp.n_term
        else:
            ge2com = ge_comp.c_term
        assert np.isclose(ge2com.value, ge_exp[2]), logging.warning(ge2com, ge_exp)


def test_ge_configuration():
    """Test the possible selections for the GE of the whole configuration."""
    ge_loops_computed = gent.ge_loops(cm_1ucs, ca_positions, 10, backend="cython")
    modes = ("max", "weighted")

    for mode in modes:
        ge_selected = gent.ge_configuration(ge_loops_computed, 0, mode)
        assert np.isclose(ge_selected.value, ge_1ucs_native_modes[mode], atol=0.0001)
        if mode == "max":
            assert ge_selected.loop[0] == ge_1ucs_native_modes["max_loop"][0]
            assert ge_selected.loop[1] == ge_1ucs_native_modes["max_loop"][1]
            assert ge_selected.thread[0] == ge_1ucs_native_modes["max_thr"][0]
            assert ge_selected.thread[1] == ge_1ucs_native_modes["max_thr"][1]
