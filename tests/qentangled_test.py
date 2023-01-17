"""
Test functions for the qentangled module
"""
import pathlib

import pytest
import numpy as np

from pyge.qentangled import qentangled

@pytest.mark.parametrize("hill_coeff, comparison", [
    (1.0, -0.1876477248),
    (1.5, -0.1812388747),
    (2.0, -0.176858283),
    (2.5, -0.174749319),
    (3.0, -0.174636300),
    (5.0, -0.184886674)
])
def test_qentangled(hill_coeff, comparison):
    """Results should be consistent with independent calculations for 1UCS"""
    parent = pathlib.Path(__file__).parent.resolve()

    ge_1ucs_native = np.load(
        parent / "data/GE_1ucs_native.npy",
        allow_pickle=True)

    # count contacts with loops and threads |j-i|>9
    n_contacts = len(
        [g for g in ge_1ucs_native if (g[0][1]-g[0][0]>9) and (g[1][1]-g[1][0]>9)]
    )

    assert n_contacts == 128

    res = qentangled(
        ge_1ucs_native,
        min_loop_len=10,
        min_thr_len=10,
        activation_params={
            "threshold": 0.7,
            "hill_coeff": hill_coeff
        }
    )

    assert np.isclose(res, comparison, atol=1e-6)
