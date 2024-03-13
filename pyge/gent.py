"""Core functions to compute the Gaussian Entanglement.

Functions to analyze topological complexes in protein configurations
through the Gaussian Entanglement as defined in:
- Baiesi et al. Sci. Rep. (2016)
- Baiesi et al. Sci. Rep. (2019)
"""

from dataclasses import dataclass
from typing import List

import numpy as np
from numba import njit

from pyge.activation import hill_fun
from pyge.libcython.cython_gaussian_entanglement import cython_gaussian_entanglement

# Global parameters
back_allowed = ("numpy", "cython")
mode_allowed = ("max", "weighted")


def _loop_list(contact_map):
    """List of loops indexes.

    Return list of two numbers, each defining the beginning and the
    end of a loop in a single polypeptide chain

    Parameters
    --------
    contact_map : array_like
        2D numpy array where entry (i,j) is !=0 if there is a contact
        between residue i and j (counting starts from 0)

    Returns
    -------
    Array : array_like
        2D array with shape (# of contacts, 2) that describe the list of
        contacting residues
    """
    # select only upper triangular matrix
    cm_triu = np.triu(contact_map, 1)
    # in a single polypeptide chain, each contact creates a loop
    idx = np.where(cm_triu != 0)
    return np.hstack((idx[0][..., None], idx[1][..., None]))


def _midposition_vectors(bead_position):
    """Compute R vectors for GE calculation, i.e. bead positions."""
    return ((bead_position[1:] + bead_position[:-1]) / 2).astype(float)


def _bond_vectors(bead_position):
    """Compute DeltaR for GE calculation, i.e. bond vectors."""
    return (bead_position[1:] - bead_position[:-1]).astype(float)


@njit
def gaussian_entanglement(loop0, loop1, thr0, thr1, positions, bonds):
    """G' one loop.

    Gaussian entanglement for one loop starting from residue i1 to i2 and
    thread starting from residue j1 to j2.

    Pure python implementation of the corresponding cython function
    See gaussian_entanglement.computeGE() for more info
    """
    # checks
    if loop0 > loop1 or thr0 > thr1:
        raise ValueError("Incorrect values for the init and end of loop or thread")

    gauss = 0
    # in baiesi2019 paper is loop1-1, but range already stops at loop1-1
    for i in range(loop0, loop1):
        dR = positions[i] - positions[thr0:thr1]
        ddR = np.cross(bonds[i], bonds[thr0:thr1])
        gauss += np.sum(
            np.sum((dR * ddR), axis=1) / np.power(np.sum(dR**2, axis=1), 3 / 2)
            # dR*ddR shape = (10, 3)
        )

        # Readable but slow for loop
        # for j in range(thr0, thr1):
        #     dR = positions[i] - positions[j]
        #     ddR = np.cross(bonds[i], bonds[j])
        #     gauss += (dR/(dR.dot(dR)**(3/2))).dot(ddR)

    return gauss / (4 * np.pi)


@dataclass
class GE:
    """Base object to save results of a GE calculation.

    Attributes
    ----------
        loop: List[int]
            list of two indexes referring to the starting
            (starting from 0 and the N terminal) and ending
            residues of the loop.
            Practically, those that are in contact.
        thread: List[int]
            As above, two indexes describing the starting and
            ending residues of the thread.
        value: float
            Gaussian entanglement value for that loop and thread.
            Corresponding to the maximum modulus among all threads
    """

    loop: List[int]
    thread: List[int]
    value: float


@dataclass
class GETermini:
    """Store GE result for the N and the C terminus."""

    n_term: GE
    c_term: GE


def ge_loops(contact_map, bead_position, thr_min_len, backend="cython"):
    """
    Compute G' for all loops, i.e. contacts, in a single polypeptide chain.

    Parameters
    ---------
        contact_map : array_like
            2D numpy array, contact map
        bead_position : array_like
            2D numpy array: list of the 3 dimensional position for
            each C_alpha (or residue in coarse-grained fashion).
            Assumed to be ordered from N to C terminus (as the contact_map).
        thr_min_len : int
            Let i and j be indexes for the thread;
            Then this argument set the condition:
                |j-i| =>  thr_min_len
        backend : str, default, 'cython'
            Specifies the backend for the GE calculation.
            Possible values: 'numpy', 'cython'.

    Returns
    -------
        result : List[GETermini]
    """
    if backend == back_allowed[0]:
        # Numpy
        ge_func = gaussian_entanglement
    elif backend == back_allowed[1]:
        # Cython
        ge_func = cython_gaussian_entanglement
    else:
        raise ValueError(f"Backend not valid, use one of the following: {back_allowed}")

    len_chain = bead_position.shape[0]
    if contact_map.shape[0] != len_chain:
        raise ValueError("Contact map and bead position have different length")
    loops = _loop_list(contact_map)
    # note: first index is 0
    # note: second index > first index, for construction

    midpos = _midposition_vectors(bead_position)  # len(midpos) = len chain - 1
    bonds = _bond_vectors(bead_position)

    result = []
    for loop in loops:
        # loop extrema
        i1 = loop[0]
        i2 = loop[1]

        GE_max_n = 0
        j1_max_n = 0
        j2_max_n = 0
        GE_max_c = 0
        j1_max_c = 0
        j2_max_c = 0
        # Search before the loop: N terminus
        for j1 in range(i1 - thr_min_len):
            for increment in range(i1 - thr_min_len - j1):
                j2 = j1 + thr_min_len + increment
                GE_loop = ge_func(i1, i2, j1, j2, midpos, bonds)
                if abs(GE_loop) > abs(GE_max_n):
                    GE_max_n = GE_loop
                    j1_max_n = j1
                    j2_max_n = j2
        # search after the loop
        for j1 in range(i2 + 1, len_chain - thr_min_len):
            for increment in range(len_chain - thr_min_len - j1):
                j2 = j1 + thr_min_len + increment
                GE_loop = ge_func(i1, i2, j1, j2, midpos, bonds)
                if abs(GE_loop) > abs(GE_max_c):
                    GE_max_c = GE_loop
                    j1_max_c = j1
                    j2_max_c = j2
        # loop's GE is the maximum in modulus
        result.append(
            GETermini(
                n_term=GE(loop=(i1, i2), thread=(j1_max_n, j2_max_n), value=GE_max_n),
                c_term=GE(loop=(i1, i2), thread=(j1_max_c, j2_max_c), value=GE_max_c),
            )
        )
    return result


def _max_between_termini(ge_termini):
    """Max between GE of the N and C terminus."""
    if abs(ge_termini.n_term.value) >= abs(ge_termini.c_term.value):
        return ge_termini.n_term
    else:
        return ge_termini.c_term


def ge_configuration(ge_list_complete, loop_min_len, mode="max", **kwards):
    r"""GE of a chain configuration with corresponding loop-thread.

    The user can chose the mode use to select the Gaussian Entanglement
    for the whole configuration.
    max:
        Return the GE that maximize the absolute value for all GE.
        Returning the corresponding loop and thread.
    weighted:
        Return a weighted average of the GE for all loops, where the
        weights are computed using an hill function of the absolute
        value of the GE. The equation is:
            \sum_{i,j \in C} G^\prime (i,j)
            \frac{f(|G^\prime| | g_0, w)}{\sum_{i,j \in C} f(|G^\prime| | g_0, w)}
        where C is the set of formed native contacts.

    Parameters
    ---------
        ge_list_complete : List[GETermini]
            Output of gent.ge_loops
        loop_min_len : int
            Let i and j be indexes for the thread;
            Then this argument set the condition:
                |j-i| =>  loop_min_len
        mode : str, default, 'max'
            how to select the G' for the whole configuration.
            The possible ones are: 'max' or 'weighted'
        kwargs:
            activation_params : dict, optional
                Parameters passed to the activation Hill function.
                threshold : float, default 0.5
                    g_0 parameter
                hill_coeff : float, default 3.0
                    Exponent of the Hill function

    Returns
    ---------
        out : GE
            GE object.
            It returns None if no loop is found satisfying the m0 threshold.
            For the weighted case, the loop and the thread are set to None
    """
    if mode == mode_allowed[0]:
        # Maximum
        fun_selection = max
        params = {"key": abs}
    elif mode == mode_allowed[1]:
        activation_function = hill_fun
        # get exponent if present, otherwise set it to 3 (see salicari2023)
        hill_coeff = kwards.get("hill_coeff", 3)
        threshold = kwards.get("threshold", 0.5)
        activation_params = {"threshold": threshold, "hill_coeff": hill_coeff}
        # there is no loop or thread associated to the average value
        loop, thread = None, None
    else:
        raise ValueError(
            f"Mode {mode} not valid, use one of the following: {mode_allowed}"
        )

    ge_list_single_terminus = [_max_between_termini(x) for x in ge_list_complete]
    ge_loop_filtered = [
        x for x in ge_list_single_terminus if x.loop[1] - x.loop[0] >= loop_min_len
    ]
    if len(ge_loop_filtered) == 0:
        return None
    ge_values = [ge.value for ge in ge_loop_filtered]

    if mode == mode_allowed[0]:
        ge_selected = fun_selection(ge_values, **params)
        loop = ge_loop_filtered[ge_values.index(ge_selected)].loop
        thread = ge_loop_filtered[ge_values.index(ge_selected)].thread
    elif mode == mode_allowed[1]:
        weights = activation_function(np.abs(ge_values), **activation_params)
        ge_selected = np.sum(weights * ge_values) / np.sum(weights)

    return GE(loop=loop, thread=thread, value=ge_selected)
