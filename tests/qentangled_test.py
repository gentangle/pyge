"""
Test functions for the qentangled module
"""
import pytest
import numpy as np

from pyge.qentangled import qentangled

# 1 riprodurre condizione dei file di antonio, ovvero modificare la lista di contatti
# per 1ucs e qualche altra proteina
# 2. calcolo con la mia funzinoe e asser

@pytest.mark.parametrize("input, thr, hill_coeff", [
    (0.0001, 0.5, 1),
    (0.0001, 0.5, 4)
])
def hill_fun_left_extrema(imp, thr, hill_coeff):
    # test extrema
    assert np.isclose(qentangled.hill_fun(imp, thr, hill_coeff), 0)
