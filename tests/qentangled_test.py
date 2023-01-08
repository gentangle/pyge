"""
Test functions for the qentangled module
"""
import pathlib

import pytest
import numpy as np

from pyge.qentangled import hill_fun

@pytest.mark.parametrize("imp, thr, hill_coeff", [
    (0.001, 0.5, 2),
    (0.001, 0.5, 4),
    (0.001, 0.7, 2),
    (0.001, 0.7, 4)
])
def test_hill_fun_left_extrema(imp, thr, hill_coeff):
    """Test extrema: result should tend to 0"""
    assert hill_fun(imp, thr, hill_coeff) < 1e-4

@pytest.mark.parametrize("imp, thr, hill_coeff", [
    (1.5, 0.5, 4),
    (2, 0.7, 4)
])
def test_hill_fun_right_extrema(imp, thr, hill_coeff):
    """Test extrema: result should tend to 1"""
    assert hill_fun(imp, thr, hill_coeff) > 0.95

def test_qentangled():
    """Results should be consistent with independent calculations"""
    parent = pathlib.Path(__file__).parent.resolve()

    ge_1ucs_native = np.load(
        parent / "data/GE_1ucs_native.npy",
        allow_pickle=True)


# 1 riprodurre condizione dei file di antonio, ovvero modificare la lista di contatti
# per 1ucs e qualche altra proteina
# 2. calcolo con la mia funzinoe e asser
