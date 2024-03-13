"""Test functions for the trajectory module."""

import json
import math
import pathlib

import numpy as np

from pyge import trajectory

parent = pathlib.Path(__file__).parent.resolve()
with open(parent / "data/1ucs_ge_traj.json", "r", encoding="utf-8") as fin:
    ge_1ucs_traj = json.load(fin)


def test_trajectory():
    """Test for the trajectory module."""
    # preparing the mask to select the configuration to test
    mask = np.zeros(41667, dtype=np.bool_)
    mask[0] = 1
    mask[20000] = 1
    mask[-1] = 1

    ge_ts, cm = trajectory.trajectory(
        parent / "data" / "1ucs_topo_CA.pdb",
        parent / "data" / "1ucs_traj.dcd",
        mask,
        {"thr_min_len": 10, "whole_config": True, "loop_min_len": 10, "mode": "max"},
        {
            "pdb_for_cm": parent / "data" / "1ucs.pdb",
            "model_id": 1,
            "chain_id": "A",
            "threshold": 4.5,
            "to_ignore": ["HOH"],
            "gamma": 1.2,
        },
    )

    assert len(ge_ts) == len(ge_1ucs_traj["results"])
    for ge, ge_res in zip(ge_ts, ge_1ucs_traj["results"]):
        assert math.isclose(ge.value, ge_res[2])
        assert ge.loop[0] == ge_res[0][0]
        assert ge.loop[1] == ge_res[0][1]
        assert ge.thread[0] == ge_res[1][0]
        assert ge.thread[1] == ge_res[1][1]

    mask = [0, 20000, -1]

    ge_ts = trajectory.trajectory(
        parent / "data" / "1ucs_topo_CA.pdb",
        parent / "data" / "1ucs_traj.dcd",
        mask,
        {"thr_min_len": 10, "whole_config": True, "loop_min_len": 10, "mode": "max"},
        {
            "model_id": 1,
            "chain_id": "A",
            "threshold": 4.5,
            "to_ignore": ["HOH"],
            "gamma": 1.2,
        },
    )

    mask = {"start": 0, "stop": 41667, "step": 20000}

    ge_ts = trajectory.trajectory(
        parent / "data" / "1ucs_topo_CA.pdb",
        parent / "data" / "1ucs_traj.dcd",
        mask,
        {"thr_min_len": 10, "whole_config": True, "loop_min_len": 10, "mode": "max"},
        {
            "model_id": 1,
            "chain_id": "A",
            "threshold": 4.5,
            "to_ignore": ["HOH"],
            "gamma": 1.2,
        },
    )


if __name__ == "__main__":
    test_trajectory()
