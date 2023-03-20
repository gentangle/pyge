"""Standard protein letter codes.

File forked and adapted from the 1.79 version of
Biopython.Bio.Data.IUPACData
"""

protein_letters_3to1 = {
    'ALA': 'A',
    'CYS': 'C',
    'ASP': 'D',
    'GLU': 'E',
    'PHE': 'F',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LYS': 'K',
    'LEU': 'L',
    'MET': 'M',
    'ASN': 'N',
    'PRO': 'P',
    'GLN': 'Q',
    'ARG': 'R',
    'SER': 'S',
    'THR': 'T',
    'VAL': 'V',
    'TRP': 'W',
    'TYR': 'Y'
}

# Here the extended AA are described
#   B = "Asx";  aspartic acid or asparagine (D or N)
#   X = "Xxx";  unknown or 'other' amino acid
#   Z = "Glx";  glutamic acid or glutamine (E or Q)
#   http://www.chem.qmul.ac.uk/iupac/AminoAcid/A2021.html#AA212
#
#   J = "Xle";  leucine or isoleucine (L or I, used in NMR)
#   Mentioned in http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html
#   Also the International Nucleotide Sequence Database Collaboration (INSDC)
#   (i.e. GenBank, EMBL, DDBJ) adopted this in 2006
#   http://www.ddbj.nig.ac.jp/insdc/icm2006-e.html
#
#   Xle (J); Leucine or Isoleucine
#   The residue abbreviations, Xle (the three-letter abbreviation) and J
#   (the one-letter abbreviation) are reserved for the case that cannot
#   experimentally distinguish leucine from isoleucine.
#
#   U = "Sec";  selenocysteine
#   http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html
#
#   O = "Pyl";  pyrrolysine
#   http://www.chem.qmul.ac.uk/iubmb/newsletter/2009.html#item35

protein_letters_3to1_extended = {
    'ALA': 'A',
    'CYS': 'C',
    'ASP': 'D',
    'GLU': 'E',
    'PHE': 'F',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LYS': 'K',
    'LEU': 'L',
    'MET': 'M',
    'ASN': 'N',
    'PRO': 'P',
    'GLN': 'Q',
    'ARG': 'R',
    'SER': 'S',
    'THR': 'T',
    'VAL': 'V',
    'TRP': 'W',
    'TYR': 'Y',
    'ASX': 'B',
    'XAA': 'X',
    'GLX': 'Z',
    'XLE': 'J',
    'SEC': 'U',
    'PYL': 'O'
}
