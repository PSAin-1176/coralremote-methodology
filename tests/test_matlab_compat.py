import numpy as np
import pytest

from coralremote.matlab_compat import find_colmajor, mquantile, mround, poly2mask


def test_mround_half_away_from_zero():
    assert list(mround([0.5, 1.5, 2.5, -0.5, -2.5, 2.4999])) == [1, 2, 3, -1, -3, 2]


def test_mquantile_matches_matlab_definition():
    x = np.array([3.0, 1.0, 2.0, 4.0])
    # MATLAB: quantile([1 2 3 4], [0.1 0.25 0.5 0.8 0.95]) = 1, 1.5, 2.5, 3.7, 4
    got = [mquantile(x, p) for p in (0.1, 0.25, 0.5, 0.8, 0.95)]
    assert np.allclose(got, [1.0, 1.5, 2.5, 3.7, 4.0])
    assert np.isnan(mquantile(np.array([np.nan]), 0.5))
    assert mquantile(np.array([1, np.nan, 3.0]), 0.5) == 2.0


@pytest.mark.parametrize("x0,y0,w,h", [(1, 1, 2, 2), (3, 2, 4, 1), (2, 5, 1, 3)])
def test_poly2mask_integer_rectangle_area(x0, y0, w, h):
    # Integer-vertex rectangles: MATLAB keeps exactly w*h pixels, the ones
    # right of / below the top-left edges.
    x = [x0, x0 + w, x0 + w, x0]
    y = [y0, y0, y0 + h, y0 + h]
    bw = poly2mask(x, y, 12, 12)
    assert bw.sum() == w * h
    rr, cc = np.nonzero(bw)
    assert rr.min() == y0 and cc.min() == x0          # 0-based = 1-based + 1 - 1


def test_poly2mask_direction_invariant():
    x = np.array([2, 7.3, 5.1, 1.2])
    y = np.array([1.4, 2.2, 8.9, 6.0])
    assert np.array_equal(poly2mask(x, y, 10, 10), poly2mask(x[::-1], y[::-1], 10, 10))


def test_poly2mask_clips_to_grid():
    bw = poly2mask([-5, 20, 20, -5], [-5, -5, 3, 3], 6, 6)
    assert bw[:3].all() and not bw[3:].any()


def test_find_colmajor_order():
    m = np.array([[1, 0, 1], [1, 1, 0]], bool)
    r, c = find_colmajor(m)
    assert list(zip(r, c)) == [(0, 0), (1, 0), (1, 1), (0, 2)]
