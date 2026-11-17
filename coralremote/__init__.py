"""CoralRemote: PlanetScope coral bleaching pipeline (Python).

Modules: ``io``, ``paths`` and ``units`` (data access and the unit contract);
``preprocess`` (Phase 3-5 rebuild for Stage 1); ``features``, ``stage1``,
``stage2``, ``train``, ``cart`` and ``forest`` (classification); ``mask``
(baseline coral mask); ``wave`` (FFT sea-state diagnostic); ``evaluate``;
``regression`` (P8 test tile). ``matlab_compat`` and ``rasterize`` hold the
MATLAB-exact helpers used for paper-of-record parity.
"""
__version__ = "1.0.0"
