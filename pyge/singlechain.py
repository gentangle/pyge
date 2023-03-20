"""GE for a single chain.

Module to contain functions used as addition layer to handle
topology files, such as PDBs, and compute Gaussian Entanglements
of single structures.

These functions do not provide tools to check the completeness and
the correctness of the structure used for the computations.
The user is responsible for this matter.
"""
import logging
import re
import sys
from dataclasses import dataclass
from typing import List

import MDAnalysis as mda

from pyge import gent
from pyge.contacts.contactmap import compute_contactmap
from pyge.gent import GE, GETermini


logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


@dataclass
class GEChain:
    """Complete GE result for a single chain."""

    loops: List[GETermini]
    global_max: GE
    global_weighted: GE


def _ca_from_topology(topology_file, selection_options):
    """CA position vectors from the PDB (topology) file.

    Parameters
    ----------
    topology_file : str or Path
        Path to the PDB file
    selection_options : str
        String to be used in the selection of atoms. For example, to select
        only the alpha carbon atoms of a specific chain, the string should be
        "and chain A" (note, the library should distinguish between
        alpha carbon and Calcium ions). See MDAnalysis documentation for more details

    Returns
    -------
    out : array_like
        Array of shape (N, 3) containing the position vectors of the alpha carbon atoms
    """
    if selection_options != "":
        selection_options = " " + selection_options
    universe = mda.Universe(str(topology_file))
    ca_positions = universe.select_atoms("name CA" + selection_options).positions
    # TODO: for now the function prints the number of CAs
    # with the aim that the user check if the selection is correct
    logging.debug(f"Number of CA atoms selected: {len(ca_positions)}")
    return ca_positions


def ge_from_pdb(pdb_file, ge_options, cm_options, selection_options=None):
    """Gaussian entanglement from a PDB file.

    This function return a GEResult object containing the list
    of Gaussian entanglements for each contact (a.k.a. loop), the maximum
    Gaussian entanglement and the weighted Gaussian entanglement.
    Starting from the PDB file.

    This function uses the ATOM coordinates from the pdb file,
    hence is responsibility of the user to check the correctness of the file.

    The user control the GE calculation through the ge_options dictionary
    and the contacts calculation through the cm_options dictionary.
    The `selection_options` is a string that is passed to the PDBParser
    to select only certain atoms. Chain and model selection are controlled
    through the cm_options dictionary.
    The `altloc` option can be modified by the user through
    the `selection_options` and `cm_options` dictionary.
    If so, the function will give priority to the `selection_options`.

    Example:
    ge_from_pdb(
        "1srl.pdb",
        ge_options={"thr_min_len":10, "loop_min_thr":0},
        cm_options={
            "model_id":1, "chain_id":"A", "threshold":4.5, "to_ignore":["HOH"]
            }
        selection_options="and altloc A"
    )

    Parameters
    ----------
    pdb_file : str or Path
        path to PDB file
    ge_options : Dict
        Dictionary containing the options for the GE calculation.
            thr_min_len : int
                Minimum length of a thread
            loop_min_len : int
                Minimum length of a loop
    cm_options : Dict
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
    selection_options : str, optional
        Selection grammar to select atoms from the PDB file.
        See https://userguide.mdanalysis.org/stable/selections.html
        for more details

    Returns
    -------
    GEChain
        Object containing the results of the GE calculation.
        The object contains the following attributes:
            loops : List[GETermini]
                List of Gaussian entanglements for each contact (a.k.a. loop)
                where the first list contains starting and finishing residue indexes
                of the loop, the second list contains starting and finishing residue
                of the thread (both starting from 0) and the float is the GE value.
            global_max : GE
                Same structure as above for the maximum GE value of the whole chain
            global_weighted : GE
                Same structure as above for the weighted GE of the whole chain
    """
    # Setup variables
    if "model_id" in cm_options:
        model_id = cm_options["model_id"]
    else:
        model_id = 1
        logging.warning(
            "WARNING: you did not provide the model_id, using the default value (1)"
        )
    if "chain_id" in cm_options:
        chain_id = cm_options["chain_id"]
    else:
        raise ValueError("You must provide the keyword 'chain_id'")
    if "threshold" in cm_options:
        threshold = cm_options["threshold"]
    else:
        raise ValueError("You must provide the keyword 'threshold'")
    if "to_include" in cm_options:
        to_include = cm_options["to_include"]
    else:
        to_include = None
    if "to_ignore" in cm_options:
        to_ignore = cm_options["to_ignore"]
    else:
        to_ignore = None

    # default altloc if not modified by the user
    altloc = "A"
    if selection_options is not None:
        if "altloc" in selection_options:
            search = re.search("altloc", selection_options)
            altloc = selection_options[search.span()[1]: search.span()[1] + 2].strip(
                " "
            )
            assert len(altloc) == 1
        elif "altloc" in cm_options:
            altloc = cm_options["altloc"]
        selection = selection_options
    else:
        selection = ""
        if "altloc" in cm_options:
            altloc = cm_options["altloc"]

    pdb_file = str(pdb_file)
    ca_positions = _ca_from_topology(pdb_file, selection)
    cm = compute_contactmap(
        pdb_file,
        model_id,
        chain_id,
        threshold,
        altloc=altloc,
        to_include=to_include,
        to_ignore=to_ignore,
    )

    if "thr_min_len" in ge_options:
        thr_min_len = ge_options["thr_min_len"]
    else:
        thr_min_len = 0
    if "loop_min_len" in ge_options:
        loop_min_len = ge_options["loop_min_len"]
    else:
        loop_min_len = 0

    ge_loops_result = gent.ge_loops(cm, ca_positions, thr_min_len)
    ge_config_max = gent.ge_configuration(ge_loops_result, loop_min_len, "max")
    ge_config_w = gent.ge_configuration(
        ge_list_complete=ge_loops_result,
        loop_min_len=loop_min_len,
        mode="weighted",
        kwards={"hill_coeff": 3, "threshold": 0.5},
    )

    return GEChain(
        loops=ge_loops_result, global_max=ge_config_max, global_weighted=ge_config_w
    )


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
#         Select the GE for the whole configuration.
#         See gent.ge_configuration for more details
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
