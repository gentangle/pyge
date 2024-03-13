"""GE for a trajectory.

Module to contain functions used as addition layer to handle
trajectory files, such as DCD files, and compute Gaussian Entanglements
timeseries.

These functions do not provide tools to check the completeness and
the correctness of the structure used for the computations.
The user is responsible for this matter.
"""

import logging
import sys

from pyge import gent
from pyge.contacts import check_formed
from pyge.singlechain import _ca_selection_from_topology

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


def trajectory(
    topology_file,
    trajectory_file,
    mask,
    ge_params,
    contacts_params,
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
    mask : array_like
        boolean 1D array with the same dimension of the timeseries used to
        select which configuration to sample
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
        - contact_map : array_like
            Used for calculation of GE native. Native contact map used to select
            the threshold for defining a contacts (see gamma)
        - gamma : float
            multiplicative factor used to select a contact threshold
            A contact is formed if distance[i,j] <= gamma*contact_map[i,j]


    Return
    ------
        ge_timeseries : List[GETermini] or List[GE]
            If ge_params['whole_config'] is True, then only
            info about the global GE is saved (List[GE]).
            Otherwise, all info are kept (List[GETermini]).
    """
    universe, ca_selection = _ca_selection_from_topology(
        topology_file, selection_options, trajectory_file
    )

    ge_timeseries = []
    for idx, _ in enumerate(universe.trajectory):
        # loop over the iterable Universe.trajectory. This is inherent from
        # the Reader class which contains ONE SNAPSHOT at a time
        if mask[idx]:
            ca_positions = ca_selection.positions

            contact_map_config = check_formed.cm_formed(
                contacts_params["contact_map"], ca_positions, contacts_params["gamma"]
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

    return ge_timeseries
