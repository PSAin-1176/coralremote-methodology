# Release checklist for v1.0.0

**Status (21 Nov 2026): prepared locally on branch `r1-release`; nothing has been pushed, tagged or published.**
- The release branch is `r1-release` (from `python-port`, plus the R1 runner, Table 3 script, configurations and rewritten README in `r1/`). `python-port` is unchanged.
- `r1-release` fast-forwards from GitHub `main` (2c0d556), so no history is rewritten on GitHub.
- The `github` remote is configured for fetching. Pushing needs the author's go-ahead.
- Still open: Codex pre-push review of the `r1-release` diff, and the author's go-ahead.

## Done
- [x] **No data tracked.** `git ls-files | grep -Ei '\.(mat|tif|tiff|npy|npz|h5|geojson|shp|pkl)$'` prints nothing. `tests/golden/regression_tile.json` holds only aggregates (400-pixel means, fractions, class counts), the same kind of numbers the paper publishes.
- [x] `paths.local.json` is untracked (`.gitignore`). No active machine-specific path defaults remain in code: path settings come from `paths.local.json` or environment variables. Explanatory references to the original `C:/SAT` locations remain in docstrings and provenance notes, and `paths.example.json` contains placeholder Windows paths. The Allen Coral Atlas path is now a setting (`aca_bathymetry`).
- [x] **Commit identity:** every commit after `main` uses the author's academic address. The three `main` commits are already on GitHub and are not rewritten.
- [x] **README:**
  - the paper's results came from the author's original MATLAB workflow (surviving scripts and archived outputs listed in `docs/PROVENANCE.md`, not in this repository); `legacy/matlab/` is an earlier public two-class workflow, not the paper's code;
  - Python is the verified reimplementation (`compat` reproduces Tables 4–5; `fixed` corrects the input scaling);
  - the units numbers match the paper's paired comparison (§3.5.6).
- [x] **Version 1.0.0:** set in `pyproject.toml` and `coralremote/__init__.py`. The CHANGELOG is updated, and `CITATION.cff` has been added.
- [x] **Internal port reports kept** (`docs/PY4–PY8`, `R1_LOO_DATE`). Local paths were removed. Their R1 details are the ones the revised paper discloses.
- [x] **Pre-push review (Codex, 18 Nov 2026) R1–R7 applied:**
  - provenance wording corrected (R1);
  - input manifest with checksums in `docs/INPUTS.md` (R2);
  - stack cache keyed on and validated against the source file (R3);
  - unit check enforced where Stage 2 reads the stack (R4);
  - inaccessible optional data skips instead of failing (R5);
  - output arguments checked with `paths.assert_writable` (R6);
  - environment and checklist wording tidied (R7).
- [x] **No second-manuscript material:** the change-detection work lives in a separate local repository, `coralremote-change`, which is not part of this release.

## R1 additions
- [x] `r1/` holds the R1 runner, `legacy_replay.py`, `table3_retest.py`, three configurations and `r1/README.md`. Configuration paths are placeholders; every SHA-256 and size is unchanged.
- [x] README results framing: R1 tables come from `r1/`; `compat` reproduces the submitted tables (provenance); `fixed` corrects the input scaling.
- [x] `docs/INPUTS.md` lists StableSand and the coregistered UDM2 flag files by hash and size.
- [ ] Codex pre-push review of the `r1-release` diff.
- [ ] Items marked unverified in `r1/README.md` checked, or accepted as unverified.

## Push (author)
1. Push only after the R1 manuscript is final and the author gives the go-ahead.
2. **Connect the repository to Zenodo first** (zenodo.org, GitHub integration, switch the repository on). A release made before this mints no DOI.
3. Push the release branch (fast-forward, no force):

   ```bash
   git push github r1-release:main
   ```

   Alternatively push `r1-release` and open a PR into `main` on GitHub.
4. Tag and release from `r1-release`, once Zenodo is connected:

   ```bash
   git tag -a v1.0.0 -m "CoralRemote 1.0.0 (AJSTR R1)"
   git push github v1.0.0
   ```

   Then create the GitHub release from the tag (Zenodo mints the DOI) and date the CHANGELOG entry.
5. Add the DOI badge to the README (one more commit and push). Add the repository URL and DOI to the paper's Data Availability statement.
