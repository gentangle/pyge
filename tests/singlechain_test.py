"""
Test functions for the singlechain module
"""
import pathlib
import json

import numpy as np

from pyge.singlechain import ge_from_pdb

parent = pathlib.Path(__file__).parent.resolve()
with open(parent / "data/1ucs_ge_configurations_modes.json", "r", encoding="utf-8") as fin:
    ge_1ucs_native_modes = json.load(fin)

# load ge_results for the two test cases
ge_1ucs_native = np.load(
    parent / "data/GE_1ucs_native.npy",
    allow_pickle=True)

def test_ge_from_pdb():
    # Test case 1: expected output provided
    expected_output_1 = {
        "loop_thr_ge": ge_1ucs_native,
        "ge_max": [
            ge_1ucs_native_modes["max"],
            ge_1ucs_native_modes["max_loop"],
            ge_1ucs_native_modes["max_thr"]
        ],
        "ge_weighted": [
            None,
            None,
            -0.68
        ]
    }
    result_1 = ge_from_pdb(
        parent / "data" / "1ucs.pdb",
        ge_options={"thr_min_len": 10, "loop_min_thr": 0},
        cm_options={
        "model_id": 1, "chain_id": "A", "threshold": 4.5, "to_ignore": ["HOH"]},
        selection_options=None
    )
    assert result_1["loop_thr_ge"] == expected_output_1["loop_thr_ge"]
    assert result_1["ge_max"] == expected_output_1["ge_max"]
    assert result_1["ge_weighted"] == expected_output_1["ge_weighted"]
