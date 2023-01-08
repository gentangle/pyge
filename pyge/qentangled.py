"""
Functions to compute the fraction of entangled contacts

The module implement the function to compute the fraction of
entangled contacts and this value from a timeseries.
Differently from trajectory and gent, here the functions takes
as input the list produced by gent.ge_loops

The function takes as argument the list of contacts together with info
about the respective thread and G' value. This list is the one
produced by gent.ge_loops
"""
import numpy as np

def hill_fun(inp, threshold, hill_coeff) -> np.ndarary:
    """
    Hill function

    Parameters
    ----------
    inp : array_like
        float or array for the x value
    threshold : float
        value characterizing the cutoff for the function
    hill_coeff : float or int
        exponent of the function

    Returns
    -------
    np.ndarray
        hill function for the parameters provided
    """
    return 1/( 1 + (threshold/inp)**hill_coeff )

def qentangled(contact_thr_ge_list, n_contacts, min_loop_len=0, activation_fun=hill_fun, activation_params=None) -> float:
    r"""
    Compute the fraction of entangled contacts defined as:
        \frac{1}{N} \sum_{formed contacts} \sgn (G') f( |G'| )

    Parameters
    ----------
    contact_thr_ge_list : see gent.ge_loops
        List of loops, relative thread and G' value
    n_contacts : int
        total number of contacts
    min_loop_len : int, optional
        Minimal number of residues to identify a chain segment as a loop, by default 0
    activation_fun : function, optional
        Python function used as activation one, by default hill_fun
    activation_params : dict, optional
        Parameters passed to the activation function;
        The input of the function must be the first argument, while
        those passed through this dict are subsequent
        By default None => {"thr":0.7, "hill_coeff":4}

    Returns
    -------
    float
        Fraction of entangled contacts
    """
    if activation_params is None:
        activation_params = {"thr":0.7, "hill_coeff":4}

    # extract ge filtering wrt loop len. Default is no filtering
    ge_array = np.array([x for x in contact_thr_ge_list if x[0][1]-x[0][0] >= min_loop_len])

    # np.sign -> 0 if the argument is 0. Ok because otherwise the
    # activation function would be 0 (and it is test to be in this way)
    # Moreover, the activation function takes as input the abs of G'
    integrand = np.sign(ge_array)*activation_fun(np.abs(ge_array), **activation_params)
    return integrand.sum()/n_contacts
