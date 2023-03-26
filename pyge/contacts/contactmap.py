"""Compute and draw the interaction matrix, or contact map from a PDB structure."""
import numpy as np

from pyge.contacts.pdb_parser import get_residues


def compute_contactmap(
    file,
    model_id,
    chain_id,
    threshold,
    altloc="A",
    to_include=None,
    to_ignore=None,
    sequence=None,
):
    """Contact map extracted from a PDB structure.

    The contact definition is purely geometric.
    Two residues are said to be in contact if at least one pair of
    non-hydrogen atoms is spatially closer than the threshold.
    The matrix entry is the distance between the alpha carbons of the
    respective residues.

    Disulfide bridges are not considered.

    Parameters
    --------
        See pdb_parser.get_residues for more details
        threshold : float
            Distance below which two non-hydrogen atoms are considered in contact
        altloc : str
            Alternative location for disordered atoms. Default: 'A'
        sequence : List[str], optional
            List of 1 letter code of the protein sequence. If provided, the
            contact map will be a square matrix with the same dimensions as
            the sequence. This parameter is use to highlight gaps in the
            sequence (i.e. missing residues).

    Returns
    -------
        contact_map : array_like
            Two dimensional numpy array with entries the alpha carbon distances
            of contacting residues. All dimensions are equal to the number of
            residues in the protein
    """
    res_list = get_residues(file, model_id, chain_id, to_include, to_ignore)
    if sequence is not None:
        n_residues = len(sequence)  # this considers missing residues too
        # fill residue list with Nones corresponding to missing residues
        for idx, one_letter_code in enumerate(sequence):
            if one_letter_code == '-':
                res_list.insert(idx, None)
            else:
                # TODO: add a check that the sequence is correct
                # Idea: the user has to provide full info for consistency,
                # otherwise the program would run also with only gap info
                continue
        if len(res_list) != len(sequence):
            raise ValueError(
                (
                    "Sequence and residue list have different lengths\n"
                    f"Read: {len(res_list)}\n Seq: {len(sequence)}"
                )
            )
    else:
        n_residues = len(res_list)
    contact_map = np.zeros((n_residues, n_residues))

    # Contacts are computed between residues which are separated by at least
    # three other residues in the polypeptide chain [noel2012]
    for row in range(0, n_residues - 4):
        if res_list[row] is None:
            continue
        for col in range(row + 4, n_residues):
            if res_list[col] is None:
                continue

            res1 = res_list[row]
            res2 = res_list[col]
            stop_res1 = False
            for atom1 in res1:
                # Contact map takes into account only heavy atoms.
                # In PDBs, Hydrogen are highlighted by a string with
                # initial char equal to H. But there are atoms such as
                # Oxygen in ARG that are identified as OH
                if atom1.get_name().startswith("H"):
                    continue
                if atom1.is_disordered():
                    atom1.disordered_select(altloc)
                for atom2 in res2:
                    if atom2.get_name().startswith("H"):
                        continue
                    if atom2.is_disordered():
                        atom2.disordered_select(altloc)

                    # the minus operator between Bio.PDB.Atom.Atom
                    # objects is overloaded to return the distance
                    # between them
                    distance = atom1 - atom2
                    if distance < threshold:
                        Ca1 = res1["CA"]
                        Ca2 = res2["CA"]
                        if Ca1.is_disordered():
                            Ca1.disordered_select(altloc)
                        if Ca2.is_disordered():
                            Ca2.disordered_select(altloc)

                        contact_map[row, col] = Ca2 - Ca1
                        stop_res1 = True
                        break
                if stop_res1:
                    break
    return contact_map


if __name__ == "__main__":
    pass
