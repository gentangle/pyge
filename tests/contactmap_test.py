import csv

import numpy as np

import pyge.contacts.contactmap as cm


with open("/home/leonardo/pyge/tests/data/1srl_seqres_cut.csv") as inf:
    reader = csv.reader(inf)
    seq = list(reader)[0]

out = cm.compute_contactmap("/home/leonardo/pyge/tests/data/1srl_w_missing_res.pdb", 1, "A", 4.5, to_ignore=["HOH"], sequence=seq)

print(out.shape)
print("tail")
print(not np.any(out[0]))  # should be True
print(not np.any(out[1]))
print(not np.any(out[2]))
print(not np.any(out[3]))
print(not np.any(out[4]))
print(not np.any(out[5]))
print(not np.any(out[6]))
print(not np.any(out[7]))
# negative check
print(not np.any(out[8]))  # should be False

# added for debugging
print("middle")
print(not np.any(out[:, 59]))  # should be False
print(not np.any(out[:, 60]))  # should be True
print(not np.any(out[:, 61]))  # should be False

# Column
print("column")
print(not np.any(out[:, 59]))  # should be False
print(not np.any(out[:, 60]))  # should be True
