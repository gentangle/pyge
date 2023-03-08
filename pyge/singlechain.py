"""
Module to contain functions used as addition layer to handle
topology files, such as PDBs, and compute Gaussian Entanglements
of single structures.

These functions do not provide tools to check the completeness and
the correctness of the structure used for the computations.
The user is responsible for this matter.
"""
import re

import MDAnalysis as mda

from pyge import gent
from pyge.contactmap.contactmap import compute_contactmap

def _ca_from_topology(topology_file, selection_options):
    """
    TODO
    Handle the topology file creating an MDAnalysis.core.Universe object

    Selects only the CA atoms (note, the library should distinguish between
    alpha carbon and Calcium ions)
    """
    # print("WARNING: pyge does not check for correctness nor completeness of the chain.\
    #     Please, make sure that the topology file provided is correct\n")
    universe = mda.Universe(str(topology_file))
    ca_positions = universe.select_atoms("name CA" + selection_options).positions
    # TODO: add a check about len
    return ca_positions


# TODO serve aggiungere dei test
def ge_from_pdb(pdb_file, mode_options, parser_options=None, selection_options=None):
    """
    TODO

    If altloc is fixed in selection_option, then the keyword in parser_option is ignored for consistency
    """

    # Setup variables
    if selection_options is None:
        selection = ''
    if "model_id" in parser_options:
        model_id = parser_options["model_id"]
    else:
        raise ValueError("You must provide the keyword 'model_id'")
    if "chain_id" in parser_options:
        chain_id = parser_options["chain_id"]
    else:
        raise ValueError("You must provide the keyword 'chain_id'")
    if "threshold" in parser_options:
        threshold = parser_options["threshold"]
    else:
        raise ValueError("You must provide the keyword 'threshold'")
    if "altloc" in selection_options:
        search = re.search("altloc", selection_options)
        altloc = selection_options[search.span()[1]:search.span()[1]+2].strip(" ")
    elif "altloc" in parser_options:
        altloc = parser_options["altloc"]
    if "to_include" in parser_options:
        to_include = parser_options["to_include"]
    else:
        to_include = None
    if "to_ignore" in parser_options:
        to_ignore = parser_options["to_ignore"]
    else:
        to_ignore = None

    ca_positions = _ca_from_topology(pdb_file, selection)
    cm = compute_contactmap(
        pdb_file,
        model_id, chain_id,
        threshold,
        altloc=altloc,
        to_include=to_include,
        to_ignore=to_ignore
    )

    if "thr_min_len" in mode_options:
        thr_min_len = mode_options["thr_min_len"]
    else:
        thr_min_len = 0
    if "loop_min_len" in mode_options:
        loop_min_len = mode_options["loop_min_len"]
    else:
        loop_min_len = 0

    ge_loops_result = gent.ge_loops(cm, ca_positions, thr_min_len)
    ge_config_max = gent.ge_configuration(ge_loops_result, loop_min_len, "max")
    ge_config_w = gent.ge_configuration(ge_loops_result, loop_min_len, "weighted", {"hill_coeff":3, "thr":0.5})

    return {"loop_thr_ge": ge_loops_result, "max_ge": ge_config_max, "ge_weighted": ge_config_w}

# def singlechain(topology_file, contact_map, mode, loop_min_len=10, thr_min_len=10):
#     """
#     Compute the Gaussian Entanglement for all loops and the whole configuration
#     of a single polypeptide chain.

#     The user provides a topology file, usually a PDB, from which the alpha carbon
#     chain (backbone) is extracted to compute the self entanglement.

#     Parameters
#     ---------
#     topology_file : str or Path
#         Path to the topology file. Usually a PDB
#     contact_map : array_like
#         2D array having entries as the distance between alpha carbon belonging to
#         residues that are defined in contact by the user
#     mode : str
#         Select the GE for the whole configuration. See gent.ge_configuration for more details
#     loop_min_len : int
#         Minimum number of residues for a loop
#     thr_min_len : int
#         Minimum number of residues for a thr

#     Returns
#     -------
#     out : dict
#         Dictionary with the result of the computation. Keywords are:
#         - 'loop_thr_ge' : format of ge_loops
#             GE for all loops in the configuration
#         - 'loop_thr_ge_MODE' : format of ge_configuration
#             GE for the whole configuration as selected with MODE parameter
#     """
#     # load topology file as Universe Object
#     ca_positions = _object_from_topology(topology_file)
#     ge_loops_result = gent.ge_loops(contact_map, ca_positions, thr_min_len)
#     ge_config_result = gent.ge_configuration(ge_loops_result, loop_min_len, mode)

#     return {"loop_thr_ge": ge_loops_result, f"loop_thr_ge_{mode}": ge_config_result}
