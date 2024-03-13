"""GE for a trajectory.

Module to contain functions used as addition layer to handle
trajectory files, such as DCD files, and compute Gaussian Entanglements
timeseries.

These functions do not provide tools to check the completeness and
the correctness of the structure used for the computations.
The user is responsible for this matter.
"""

import numpy as np
from loguru import logger

from pyge import gent
from pyge.contacts import check_formed
from pyge.contacts.contactmap import compute_contactmap
from pyge.singlechain import (
    _ca_selection_from_topology,
    _check_altloc,
    _check_cm_options,
)


def trajectory(
    topology_file,
    trajectory_file,
    mask,
    ge_params,
    cm_options,
    selection_options=None,
):
    """GE for all loops or the whole configurations of single chain trajectory.

    The use case is to provide a topology and a trajectory file output of a simulation
    of a coarse-grained polypeptide chain.

    Parameters
    ---------
    topology_file : str or Path
        path to the topology file
    trajectory_file : str or Path
        path to the trajectory file
    mask : list, np.ndarray or dict
        mask to select the frames to consider. It can be a list, a numpy array
        or a dictionary with the keys: start (included), stop (excluded), step
    ge_params : dict
        Dictionary composed as follows:
        - thr_min_len : int
            minimum number of residues to have a thread
        - whole_config : bool
            If True, save *only* the Gaussian Entanglement of the whole configuration
            based on the mode keyword
        - loop_min_len : int
            minimum number of residues to have a loop
        - mode : str
            mode selection, for more details see pyge.gent.ge_configuration
    contacts_params : dict
        Dictionary containing the options for the contact map calculation.
        Each contact defines a loop
            model_id : int, optional
                Model ID for the chain (default 1, used,
                for example, for X-Ray structures)
            chain_id : str, by default 'A'
                Chain ID to select chain (auth_asym parameter in PDB)
            threshold : float
                Minimum distance for which two non-Hydrogen atoms are considered
                in contact
            to_include : List[str], optional
                List of residue names to include when parsing the PDBge_config_max
            to_ignore : List[str], optional
                List of residue names to ignore when parsing the PDB.
                E.g. water (HOH), ions etc.
            gamma : float
                multiplicative factor used to select a contact threshold
                A contact is formed if distance[i,j] <= gamma*contact_map[i,j]
            pdb_for_cm : str, optional
                Path to the PDB file to use for the contact map calculation. When specified, this is used instead
                of the topology file.
                Usually, this is used when the topo file has be
                modified (e.g. for coarse-grained models) but the
                contact map has to be calculated on all-atom model.

    Return
    ------
        ge_timeseries : List[GETermini] or List[GE]
            If ge_params['whole_config'] is True, then only
            info about the global GE is saved (List[GE]).
            Otherwise, all info are kept (List[GETermini]).
    """
    selection, altloc = _check_altloc(selection_options, cm_options)
    _check_cm_options(cm_options)

    universe, ca_selection = _ca_selection_from_topology(
        topology_file, selection, trajectory_file
    )
    if cm_options["pdb_for_cm"] is not None:
        topo_for_cm = cm_options["pdb_for_cm"]
    else:
        topo_for_cm = topology_file
    cm = compute_contactmap(
        topo_for_cm,
        cm_options["model_id"],
        cm_options["chain_id"],
        cm_options["threshold"],
        altloc=altloc,
        to_include=cm_options["to_include"],
        to_ignore=cm_options["to_ignore"],
    )

    if isinstance(mask, list):
        # case in which the mask is a list of frame to consider
        traj_slice = universe.trajectory[mask]
    elif isinstance(mask, np.ndarray):
        if mask.dtype == np.bool_:
            # case in which the mask is a boolean array selecting the frames to consider
            traj_slice = universe.trajectory[mask]
        else:
            raise ValueError("mask must be a boolean array: dtype = np.bool_")
    elif isinstance(mask, dict):
        # case in which the mask is a dictionary with teh ranges to use
        for key in mask:
            if not isinstance(mask[key], int):
                raise ValueError("mask values (start, stop, step) must be integers")
        traj_slice = universe.trajectory[mask["start"] : mask["stop"] : mask["step"]]
    else:
        raise ValueError(
            "mask must be a list, a numpy array or a dictionary with the keys: start, stop, step"
        )

    logger.info(f"Starting computation of GE for {len(traj_slice):d} frames")

    ge_timeseries = []
    for _ in traj_slice:
        ca_positions = ca_selection.positions

        contact_map_config = check_formed.cm_formed(
            cm, ca_positions, cm_options["gamma"]
        )
        ge_loops_result = gent.ge_loops(
            contact_map_config, ca_positions, ge_params["thr_min_len"]
        )

        if ge_params["whole_config"]:
            ge_config_result = gent.ge_configuration(
                ge_loops_result, ge_params["loop_min_len"], ge_params["mode"]
            )

            # save only the GE for the whole configuration
            ge_timeseries.append(ge_config_result)
            continue

        # save all GE from all loops
        ge_timeseries.append(ge_loops_result)

    return ge_timeseries, cm
