"""Compute and draw the interaction matrix, or contact map from a PDB structure."""
import numpy as np

from pyge.contacts.pdb_parser import get_residues
from pyge.contacts.protein_letters import (
    protein_letters_3to1,
    protein_letters_3to1_extended,
    )


aa_3to1 = {}
aa_3to1.update(protein_letters_3to1)
aa_3to1.update(protein_letters_3to1_extended)


def compute_contactmap(
    file,
    model_id,
    chain_id,
    threshold,
    altloc=None,
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
        altloc : str, by default: None
            Alternative location for disordered atoms. If set to None,
            the position selected are those with highest occupancy.
            See: https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
            Section: Disordered atom positions
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
            if one_letter_code == "-":
                res_list.insert(idx, None)
            else:
                # design choice: the sequence provided has to be
                # consistent with atoms entries.
                cur_code = aa_3to1[res_list[idx].get_resname()]
                if not cur_code == one_letter_code:
                    raise ValueError(
                        (
                            "Sequence and ATOM entries do not match\n"
                            f"From sequence: {one_letter_code} at idx: {idx}\n"
                            f"From ATOM: {cur_code} "
                            f"from res (id): {res_list[idx].get_id()}\n"
                            f"Full res list up to {idx} itaration: {res_list}\n"
                        )
                    )
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
                if atom1.is_disordered() and altloc is not None:
                    atom1.disordered_select(altloc)
                for atom2 in res2:
                    if atom2.get_name().startswith("H"):
                        continue
                    if atom2.is_disordered() and altloc is not None:
                        atom2.disordered_select(altloc)

                    # the minus operator between Bio.PDB.Atom.Atom
                    # objects is overloaded to return the distance
                    # between them
                    distance = atom1 - atom2
                    if distance < threshold:
                        Ca1 = res1["CA"]
                        Ca2 = res2["CA"]
                        if Ca1.is_disordered() and altloc is not None:
                            Ca1.disordered_select(altloc)
                        if Ca2.is_disordered() and altloc is not None:
                            Ca2.disordered_select(altloc)

                        contact_map[row, col] = Ca2 - Ca1
                        stop_res1 = True
                        break
                if stop_res1:
                    break
    return contact_map


if __name__ == "__main__":
    pass
