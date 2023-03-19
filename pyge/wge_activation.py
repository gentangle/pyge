"""
Module to compute [W]eighted average of [GE]
where weights are given by an Hill [A]ctivation function

--> wGEa
"""
import numpy as np

from pyge.activation import hill_fun

def wgea(contact_thr_ge_list, min_loop_len=0, min_thr_len=0, activation_fun=None, activation_params=None) -> float:
    r"""
    Compute the weighted average of G' values, where weights
    are given by the activation function passed as argument
    evaluated on |G'|.
    Definition:
        \sum_{i,j \in C} G^\prime (i,j) \frac{f(|G^\prime| | G_0, w)}{\sum_{i,j \in C} f(|G^\prime| | G_0, w)}
    where C is the set of formed native contacts

    The goal is to weight more higher values of GE

    Parameters
    ----------
    contact_thr_ge_list : see gent.ge_loops
        List of loops, relative thread and G' value
    min_loop_len : int, optional
        Minimal number of residues to identify a chain segment as a loop
        loop_end - loop_start >= min_loop_len
        by default 0
    min_thr_len : int, optional
        Minimal number of residues to identify a chain segment as a thread
        thr_end - thr_start >= min_thr_len
        by default 0
    activation_fun : function, optional
        Python function used as activation one, by default hill_fun
    activation_params : dict, optional
        Parameters passed to the activation function;
        The input of the function must be the first argument, while
        those passed through this dict are subsequent
        By default None => {"thr":0.7, "hill_coeff": 4.0}

    Returns
    -------
    float
        Weighted GE average
    """
    if activation_fun is None and activation_params is None:
        activation_fun = hill_fun
        activation_params = {"threshold": 0.5, "hill_coeff": 3}

    # extract ge filtering wrt loop len. Default is no filtering
    ge_array = np.array(
        [x[2] for x in contact_thr_ge_list if (x[0][1]-x[0][0] >= min_loop_len) and (x[1][1]-x[1][0] >= min_thr_len)]
    )
    # weights
    weights = activation_fun(np.abs(ge_array), **activation_params)

    return np.sum(weights*ge_array)/np.sum(weights)
