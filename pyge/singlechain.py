"""
Module to contain functions used as addition layer to handle
topology files, such as PDBs, and compute Gaussian Entanglements
of single structures.

These functions do not provide tools to check the completeness and
the correctness of the structure used for the computations.
The user is responsible for this matter.
"""
import MDAnalysis as mda

from pyge import gent

def _object_from_topology(topology_file):
    """
    Handle the topology file creating an MDAnalysis.core.Universe object

    Selects only the CA atoms (note, the library should distinguish between
    alpha carbon and Calcium ions)
    """
    print("WARNING: pyge does not check for correctness nor completeness of the chain.\
        Please, make sure that the topology file provided is correct\n")
    universe = mda.Universe(str(topology_file))
    return universe.select_atoms("name CA").positions


def singlechain(topology_file, contact_map, mode, loop_min_len=10, thr_min_len=10):
    """
    Compute the Gaussian Entanglement for all loops and the whole configuration
    of a single polypeptide chain.

    The user provides a topology file, usually a PDB, from which the alpha carbon
    chain (backbone) is extracted to compute the self entanglement.

    Parameters
    ---------
    topology_file : str or Path
        Path to the topology file. Usually a PDB
    contact_map : array_like
        2D array having entries as the distance between alpha carbon belonging to
        residues that are defined in contact by the user
    mode : str
        Select the GE for the whole configuration. See gent.ge_configuration for more details
    loop_min_len : int
        Minimum number of residues for a loop
    thr_min_len : int
        Minimum number of residues for a thr

    Returns
    -------
    out : dict
        Dictionary with the result of the computation. Keywords are:
        - 'loop_thr_ge' : format of ge_loops
            GE for all loops in the configuration
        - 'loop_thr_ge_MODE' : format of ge_configuration
            GE for the whole configuration as selected with MODE parameter
    """
    # load topology file as Universe Object
    ca_positions = _object_from_topology(topology_file)
    ge_loops_result = gent.ge_loops(contact_map, ca_positions, thr_min_len)
    ge_config_result = gent.ge_configuration(ge_loops_result, loop_min_len, mode)

    return {"loop_thr_ge": ge_loops_result, f"loop_thr_ge_{mode}": ge_config_result}
