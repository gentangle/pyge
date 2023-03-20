"""Module with functions used to check formed contacts during a trajectory."""
import numpy as np
from MDAnalysis.analysis import distances


# TODO: optimize this computation: I think is unnecessary to loop over the entire matrix
# because more then half of it is empty


def cm_formed(contact_map_native, ca_positions, gamma):
    """Contact map from native one and C alpha positions.

    Compute a contact map for the configuration stored in ca_positions.
    This function takes as a reference the native contact map and returns
    a boolean array in which an entry is True (hence, the contact is formed
    between two residues interaction in the native conformation) if:
        distance[i,j] <= gamma*contact_map_native[i,j]
    False otherwise.

    Parameters
    --------
        contact_map_native : array_like
            2D numpy array storing the distances between interacting residues in the
            native state. In this case, distances are between alpha carbon and the
            interaction network is chosen by the user.
            Note: only the upper triangular part is considered for the calculation
        ca_positions : array_like
            Array with shape (N,3), where N is the number of residues in the chain,
            containing alpha carbon positions
        gamma : float
            Multiplicative factor for the contact threshold (see above)

    Returns
    -------
        out : array_like
            2D array with the same dimensions as contact_map_native where each entry
            is a 1 when a *native* contact is formed, 0 otherwise
    """
    # Ensure to collect only the upper triangular matrix
    cm_native_mod = np.triu(contact_map_native)
    cm_native_mod[cm_native_mod == 0] = -1
    n_residues = contact_map_native.shape[0]

    # Compute distance matrix
    self_distances = distances.self_distance_array(ca_positions)
    sq_dist_arr = np.zeros((n_residues, n_residues))
    triu = np.triu_indices_from(sq_dist_arr, k=1)
    sq_dist_arr[triu] = self_distances

    return sq_dist_arr <= gamma * cm_native_mod
