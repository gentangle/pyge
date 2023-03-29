"""Parser to load residue list from PDB."""
import logging
import sys
from typing import List

import Bio.PDB as pdb

from pyge.contacts.protein_letters import protein_letters_3to1


logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


def get_residues(
    file, model_id, chain_id, to_include=None, to_ignore=None, debug=False
) -> List[pdb.Residue.Residue]:
    """List of Residue objects read from a PDB file.

    This function parse only standard amino acids as specified in the
    protein_letters.py. If the user does not specifies exceptions
    to this rule, a RuntimeError is raised.
    To modify this behavior, the user has to provide a list of
    strings with the ResidueIDs to be included or ignored.
    The 3 letter code is used, please provide it upper case.

    E.g.
    to_include = ['XLE', 'PYL']
    to_exclude = ['SO4', 'HOH']
    to_exclude = ['HETATM']

    The code specifies every time the parser read a user-included or
    user-excluded residue, specifying its position along the chain.

    The parser is the one implemented in the Biopython.PDB module,
    and it will raise a warning if some atoms are missing. It is up to the
    user to provide a complete PDB structure.

    Parameters
    ---------
        file : str
            Absolute/relative path to the PDB file
        model_id : int
            Selects the model (must be => 1). The same convention as in the PDB is used
        chain_id : str
            Selects the chain using the same convention as in the PDB
        to_include : List[str]
            List with the 3 letter code of residuesID to consider in the parsing apart
            from the proteinogenic ones
        to_ignore : List[str], str
            As to_include but to ignore the residue (3 letter) ID specified
            If you want to ignore water molecules, you can also provide 'W'.
            If you want to ignore all HETATM
            (hence water too),
            you can also provide 'HETATM'.

    Returns
    -------
        res_list : List[Bio.PDB.Residue.Residue]
            List of residue object
    """
    if to_include is None:
        to_include = []
    if to_ignore is None:
        to_ignore = []
    if "W" in to_ignore:
        to_ignore.append("HOH")
    if "HETATM" in to_ignore:
        ignore_hetatm_flag = True
        to_ignore.append("HOH")
    else:
        ignore_hetatm_flag = False

    parser = pdb.PDBParser()
    name_protein = file[-8:-4]  # get unique protein ID of 4 characters
    # This assumes that the file name is the protein ID
    structure = parser.get_structure(name_protein, file)

    model = list(structure.get_models())[model_id - 1]
    # -1 because Python starts indexing from 0

    chains = model.get_chains()
    found = False
    for cha in chains:
        if cha.id == chain_id:
            chain = cha
            found = True
            break
    if not found:
        raise ValueError(f"The chain ID provided {chain_id} is not valid")

    res_list = []
    for residue in chain:
        if (
            ignore_hetatm_flag
            and residue.get_id()[0].startswith("H_")
            and (residue.get_resname() not in to_ignore)
        ):
            # The last condition ensures that if the user specifies a residue,
            # e.g. a modified residue (MODRES), to be included then this condition
            # do not throw it away if HETATM is selected to be ignored.
            # Use case: PDB with ligands and modified residues, the user wants to
            # ignore ONLY the ligands.
            continue

        # filtering all proteinogenic AA
        if residue.get_resname() in protein_letters_3to1:
            res_list.append(residue)

        elif residue.get_resname() in to_include:
            if debug:
                logging.debug(
                    (
                        f"Residue {residue.get_resname()} with resID {residue.id[1]} "
                        "has been INCLUDED because specified "
                        "by the user through `to_include`"
                    )
                )
            res_list.append(residue)

        elif residue.get_resname() in to_ignore:
            if debug:
                logging.debug(
                    (
                        f"Residue {residue.get_resname()} with resID {residue.id[1]} "
                        "has been EXCLUDED because specified "
                        "by the user through `to_ignore`"
                    )
                )
            continue

        else:
            raise RuntimeError(
                (
                    f"Residue {residue.get_resname()} with resID {residue.id[1]} "
                    f"is not included in:\n{list(protein_letters_3to1.keys())}"
                )
            )
    return res_list


if __name__ == "__main__":
    pass
