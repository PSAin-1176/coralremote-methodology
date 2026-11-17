"""Py5: baseline-mask morphology chain and wave FFT diagnostic."""
import numpy as np
import pytest

from coralremote import mask, paths, wave


# ---- primitives (synthetic, no data needed) ---------------------------------
def test_disk_strel_radius1_is_cross():
    assert np.array_equal(mask.disk_strel(1), [[0, 1, 0], [1, 1, 1], [0, 1, 0]])


def test_imerode_pads_with_foreground():
    a = np.ones((4, 4), bool)
    assert mask.imerode(a, mask.disk_strel(1)).all()


def test_bwareaopen_is_8_connected():
    a = np.zeros((5, 5), bool)
    a[1, 1] = a[2, 2] = True                       # diagonal pair = one 8-connected object of 2 px
    assert mask.bwareaopen(a, 2).sum() == 2
    assert mask.bwareaopen(a, 3).sum() == 0


def test_fill_holes_4_connected_background():
    a = np.zeros((5, 5), bool)
    a[1:4, 1:4] = True
    a[2, 2] = False
    assert mask.imfill_holes(a)[2, 2]
    b = np.ones((4, 4), bool)
    b[1, 1] = b[2, 2] = False                      # diagonal background gap: not a 4-connected leak
    b[0, 0] = False
    assert mask.imfill_holes(b)[1, 1] and not mask.imfill_holes(b)[0, 0]


def test_majority_threshold_five():
    a = np.zeros((3, 3), bool)
    a.flat[[0, 1, 2, 3]] = True
    assert not mask.majority(a)[1, 1]
    a.flat[4] = True
    assert mask.majority(a)[1, 1]


def test_hann_matches_matlab_symmetric():
    h = wave._hann(5)
    assert np.allclose(h, [0, 0.5, 1, 0.5, 0])


def test_wave_annulus_bins():
    g = wave.tile_geometry()
    assert (g["rb_min"], g["rb_max"]) == (4, 64)
    assert sum(m.sum() for m in g["sector_masks"]) == g["ann"].sum()


# ---- parity on the paper-of-record data -------------------------------------
@pytest.mark.data
@pytest.mark.skipif(not (paths.available(paths.BASELINE_MASK) and paths.available(mask.PHASE9_DIR)),
                    reason="baseline mask or A3 Phase 9 maps not available")
def test_baseline_mask_chain_exact():
    import scipy.io as sio
    m = sio.loadmat(str(paths.BASELINE_MASK), squeeze_me=True)
    r = mask.build_baseline_mask(mask.load_phase9_class_maps(), m["reef_mask"] > 0)
    assert np.array_equal(r["vote_map"], m["vote_map"]) and np.array_equal(r["tier_map"], m["tier_map"])
    for k, n in (("mask0", 21460), ("maskClean", 21961), ("maskPeel", 20273), ("maskFinal", 20330)):
        assert r[k].sum() == n and np.array_equal(r[k], m[k] > 0), k
    assert np.array_equal(r["maskFinal"], m["baseline_coral_mask"] > 0)


@pytest.mark.data
@pytest.mark.skipif(not paths.available(paths.RAW_STACK), reason="raw stack not available")
def test_wave_verdicts_match_paper():
    from coralremote.preprocess import raw_stack
    s = raw_stack()
    out, _ = wave.run_wave_diagnostic(s.get, s.dates)
    assert {d: r["verdict"] for d, r in out.items()} == wave.PAPER_VERDICTS
    assert out["20230629"]["n_tiles_used"] == 38 and out["20230629"]["confidence"] == "LOW"
    assert sum(r["confidence"] == "HIGH" for r in out.values()) == 13
    # 20240514: directional AI evidence but no roughness support -> CALM (note section 7)
    assert out["20240514"]["n_bands_flag"] >= 6 and out["20240514"]["n_rough_bands"] < wave.ROUGH_MIN_BANDS
