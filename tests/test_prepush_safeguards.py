"""Data-free tests for the pre-push safeguards: cache identity (R3), the
load-time unit check (R4), optional-data detection (R5) and output-path
protection (R6)."""
import os
from pathlib import Path

import numpy as np
import pytest

from coralremote import io, paths
from coralremote.stage2 import StackView
from coralremote.units import COMPAT, FIXED, SR_SCALE, UnitError


# ---------------------------------------------------------------- R3 cache identity
def test_cache_name_depends_on_source_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("CORALREMOTE_CACHE", str(tmp_path / "cache"))
    a = tmp_path / "a" / "Stack.mat"
    b = tmp_path / "b" / "Stack.mat"
    assert io.cache_path_for(a, "stack_glint") != io.cache_path_for(b, "stack_glint")
    assert io.cache_path_for(a, "stack_glint") == io.cache_path_for(a, "stack_glint")


def test_cache_invalidated_when_source_changes(tmp_path):
    src = tmp_path / "Stack.mat"
    src.write_bytes(b"x" * 100)
    cache = tmp_path / "Stack__v.npy"
    np.save(cache, np.zeros((2, 3, 4, 8), np.float32))
    meta = io.cache_meta(src, "v", (2, 3, 4, 8), ["20230404", "20230509"])
    io.write_cache_meta(cache, meta)
    assert io.cache_is_valid(cache, meta)
    # different dates or shape than the source declares -> rebuild
    assert not io.cache_is_valid(cache, io.cache_meta(src, "v", (2, 3, 4, 8), ["20230404", "20240307"]))
    assert not io.cache_is_valid(cache, io.cache_meta(src, "v", (3, 3, 4, 8), ["a", "b", "c"]))
    # edited source (size and mtime change) -> rebuild
    src.write_bytes(b"y" * 120)
    os.utime(src, ns=(meta["mtime_ns"] + 10**9, meta["mtime_ns"] + 10**9))
    assert not io.cache_is_valid(cache, io.cache_meta(src, "v", (2, 3, 4, 8), ["20230404", "20230509"]))
    # missing sidecar -> rebuild
    io._meta_path(cache).unlink()
    assert not io.cache_is_valid(cache, meta)


# ---------------------------------------------------------------- R4 unit check at load
class _FakeStack:
    def __init__(self, values, on_disk_scale):
        self._v = values
        self.on_disk_scale = on_disk_scale
        self.dates = ["d1"]
        self.ref = None
        self.path = Path("fake.mat")

    def get(self, date):
        return self._v


def _reflectance(h=6, w=6):
    rng = np.random.default_rng(0)
    A = rng.uniform(0.02, 0.08, (h, w, 8)).astype(np.float32)
    A[..., 7] = 0.01                                   # NIR
    return A


def test_fixed_refuses_stack_already_in_reflectance():
    with pytest.raises(UnitError):
        StackView(_FakeStack(_reflectance(), on_disk_scale=1.0), FIXED)


def test_double_scaled_stack_passes_in_both_modes():
    m = np.ones((6, 6), bool)
    st = _FakeStack(_reflectance() / SR_SCALE, on_disk_scale=SR_SCALE)
    assert np.allclose(StackView(st, FIXED, check_mask=m).get("d1"), _reflectance(), rtol=1e-5)
    assert np.allclose(StackView(st, COMPAT, check_mask=m).get("d1"), _reflectance() / SR_SCALE)


def test_wrongly_declared_scale_is_caught_on_read():
    m = np.ones((6, 6), bool)
    st = _FakeStack(_reflectance(), on_disk_scale=SR_SCALE)   # claims double-scaled, is not
    with pytest.raises(UnitError):
        StackView(st, FIXED, check_mask=m).get("d1")
    with pytest.raises(UnitError):
        StackView(st, COMPAT, check_mask=m).get("d1")


# ---------------------------------------------------------------- R5 optional data
def test_available_treats_access_errors_as_missing(monkeypatch, tmp_path):
    assert paths.available(tmp_path)
    assert not paths.available(tmp_path / "absent")

    def boom(self):
        raise PermissionError("access denied")
    monkeypatch.setattr(Path, "exists", boom)
    assert not paths.available(tmp_path)


# ---------------------------------------------------------------- R6 output paths
def test_outputs_inside_protected_folders_are_refused():
    for d in paths.READ_ONLY_DIRS:
        with pytest.raises(PermissionError):
            paths.assert_writable(Path(d) / "out.json")
