"""Set of activation functions."""

import numpy as np


def hill_fun(inp, threshold: float, hill_coeff: float) -> np.ndarray:
    """
    Hill function.

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
    # use numpy array even if a single float is provided for a correct
    # handle of 0 (because we divide for 0 even thought the function tends to 0
    # for inp going to 0)
    inp = np.asarray(inp)
    with np.errstate(divide="ignore"):
        out = 1 / (1 + (threshold / inp) ** hill_coeff)
    # check fo NaN
    if np.isnan(np.sum(out)):
        raise RuntimeError("NaNs are detected as output of the Hill function!")
    return out
