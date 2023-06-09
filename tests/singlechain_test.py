"""Test functions for the singlechain module."""
import json
import pathlib

import numpy as np

from pyge.singlechain import ge_from_pdb


parent = pathlib.Path(__file__).parent.resolve()
with open(
    parent / "data/1ucs_ge_configurations_modes.json", "r", encoding="utf-8"
) as fin:
    ge_1ucs_native_modes = json.load(fin)

# load ge_results for the two test cases
ge_1ucs_native = np.load(parent / "data/GE_1ucs_native.npy", allow_pickle=True)


def test_ge_from_pdb():
    """Test the ge_from_pdb function."""
    # Test case 1: expected output provided
    expected_output_1 = {
        "loop_thr_ge": ge_1ucs_native,
        "ge_max": [
            tuple(ge_1ucs_native_modes["max_loop"]),
            tuple(ge_1ucs_native_modes["max_thr"]),
            ge_1ucs_native_modes["max"],
        ],
        "ge_weighted": [None, None, -0.68443],
    }
    result_1 = ge_from_pdb(
        parent / "data" / "1ucs.pdb",
        ge_options={"thr_min_len": 10, "loop_min_thr": 0},
        cm_options={
            "model_id": 1,
            "chain_id": "A",
            "threshold": 4.5,
            "to_ignore": ["HOH"],
        },
        selection_options=None,
    )
    for res, exp in zip(result_1.loops, ge_1ucs_native):
        if abs(res.n_term.value) >= abs(res.c_term.value):
            ge2com = res.n_term
        else:
            ge2com = res.c_term
        assert ge2com.loop == exp[0]
        assert ge2com.thread == exp[1]
        assert np.isclose(ge2com.value, exp[2])
    assert result_1.global_max.loop == expected_output_1["ge_max"][0]
    assert result_1.global_max.thread == expected_output_1["ge_max"][1]
    assert np.isclose(
        result_1.global_max.value, expected_output_1["ge_max"][2], atol=0.0001
    )
    assert np.isclose(
        result_1.global_weighted.value, expected_output_1["ge_weighted"][2], atol=0.0001
    )
