"""Test functions for the singlechain module."""

import json
import pathlib

import numpy as np
import pytest

from pyge.singlechain import GEChain, _ca_selection_from_topology, ge_from_pdb

parent = pathlib.Path(__file__).parent.resolve()
with open(
    parent / "data/1ucs_ge_configurations_modes.json", "r", encoding="utf-8"
) as fin:
    ge_1ucs_native_modes = json.load(fin)

# load ge_results for the two test cases
ge_1ucs_native = np.load(parent / "data/GE_1ucs_native.npy", allow_pickle=True)


@pytest.fixture
def setup_data_ca_sel():
    topology_file = parent / "data" / "1ucs_topo_CA.pdb"
    trajectory_file = parent / "data" / "1ucs_traj.dcd"
    return topology_file, trajectory_file


def test_ca_selection_from_topology(setup_data_ca_sel):
    topology_file, trajectory_file = setup_data_ca_sel
    universe, ca_selection = _ca_selection_from_topology(
        topology_file, "and chainID A", trajectory_file
    )
    assert universe is not None
    assert ca_selection is not None


def test_ca_selection_from_topology_no_trajectory(setup_data_ca_sel):
    topology_file, _ = setup_data_ca_sel
    universe, ca_selection = _ca_selection_from_topology(topology_file, "and chainID A")
    assert universe is not None
    assert ca_selection is not None


def test_ca_selection_from_topology_no_selection(setup_data_ca_sel):
    topology_file, _ = setup_data_ca_sel
    universe, ca_selection = _ca_selection_from_topology(topology_file)
    assert universe is not None
    assert ca_selection is not None


@pytest.fixture
def setup_data():
    # Replace with actual test data
    dummy_pdb_file = parent / "data" / "1ucs.pdb"
    dummy_ge_options = {"thr_min_len": 10, "loop_min_len": 0}
    dummy_cm_options = {
        "model_id": 1,
        "chain_id": "A",
        "threshold": 4.5,
        "to_ignore": ["HOH"],
    }
    selection_options = "and not altloc B"
    return dummy_pdb_file, dummy_ge_options, dummy_cm_options, selection_options


def test_ge_from_pdb_with_all_options(setup_data):
    pdb_file, ge_options, cm_options, selection_options = setup_data
    result = ge_from_pdb(pdb_file, ge_options, cm_options, selection_options)
    assert isinstance(result, GEChain)


def test_ge_from_pdb_without_selection_options(setup_data):
    pdb_file, ge_options, cm_options, _ = setup_data
    result = ge_from_pdb(pdb_file, ge_options, cm_options)
    assert isinstance(result, GEChain)


def test_ge_from_pdb_without_cm_options_model_id(setup_data):
    pdb_file, ge_options, cm_options, selection_options = setup_data
    cm_options.pop("model_id")
    result = ge_from_pdb(pdb_file, ge_options, cm_options, selection_options)
    assert isinstance(result, GEChain)


def test_ge_from_pdb_without_cm_options_chain_id(setup_data):
    pdb_file, ge_options, cm_options, selection_options = setup_data
    cm_options.pop("chain_id")
    with pytest.raises(ValueError):
        ge_from_pdb(pdb_file, ge_options, cm_options, selection_options)


def test_ge_from_pdb_without_cm_options_threshold(setup_data):
    pdb_file, ge_options, cm_options, selection_options = setup_data
    cm_options.pop("threshold")
    with pytest.raises(ValueError):
        ge_from_pdb(pdb_file, ge_options, cm_options, selection_options)


def test_ge_from_pdb_with_cm_options_to_include(setup_data):
    pdb_file, ge_options, cm_options, selection_options = setup_data
    cm_options["to_include"] = ["and name CA"]
    result = ge_from_pdb(pdb_file, ge_options, cm_options, selection_options)
    assert isinstance(result, GEChain)


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
