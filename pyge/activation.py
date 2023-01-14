"""
Set of activation functions
"""
import numpy as np
from numba import njit

@njit
def hill_fun(inp: np.ndarray, threshold: float, hill_coeff: float) -> np.ndarray:
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
