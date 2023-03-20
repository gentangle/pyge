"""Test function for the activation.py module."""
import pytest

from pyge.activation import hill_fun


@pytest.mark.parametrize(
    "imp, thr, hill_coeff",
    [(0.001, 0.5, 2), (0.001, 0.5, 4), (0.001, 0.7, 2), (0.001, 0.7, 4)],
)
def test_hill_fun_left_extrema(imp, thr, hill_coeff):
    """Test extrema: result should tend to 0."""
    assert hill_fun(imp, thr, hill_coeff) < 1e-4


@pytest.mark.parametrize("imp, thr, hill_coeff", [(1.5, 0.5, 4), (2, 0.7, 4)])
def test_hill_fun_right_extrema(imp, thr, hill_coeff):
    """Test extrema: result should tend to 1."""
    assert hill_fun(imp, thr, hill_coeff) > 0.95
