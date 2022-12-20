import pathlib

import numpy as np

from pyge import trajectory

parent = pathlib.Path(__file__).parent.resolve()
contact_map_native = np.load(parent / "data" / "1ucs_cm.npy")

# TODO: complete with an example run from the other code

def test_trajectory():

    mask = np.zeros(41667)
    mask[0] = 1
    mask[20000] = 1
    mask[-1] = 1
    ge_ts = trajectory.trajectory(
        parent / "data" / "1ucs_topo.pdb",
        parent / "data" / "1ucs_traj.dcd",
        mask,
        {"thr_min_len": 10, "whole_config": False},
        {"contact_map": contact_map_native ,"gamma": 1.2})

