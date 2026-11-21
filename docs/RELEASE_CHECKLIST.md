# Release checklist for v1.0.0

**Status (21 Nov 2026): prepared locally on branch `r1-release`; nothing has been pushed or published. The only tag is the local annotated `matlab-prototype` on 2c0d556 (not pushed).**
- The release branch is `r1-release` (from `python-port`, plus the R1 runner, Table 3 script, configurations and rewritten README in `r1/`). `python-port` is unchanged.
- `r1-release` fast-forwards from GitHub `main` (2c0d556), so no history is rewritten on GitHub.
- The `github` remote is configured for fetching. Pushing needs the author's go-ahead.
- Codex pre-push review of `r1-release` (21 Nov 2026) done; its seven findings are applied (see below). Still open: the author's go-ahead.

## Done
- [x] **No data tracked.** `git ls-files | grep -Ei '\.(mat|tif|tiff|npy|npz|h5|geojson|shp|pkl)$'` prints nothing. `tests/golden/regression_tile.json` holds only aggregates (400-pixel means, fractions, class counts), the same kind of numbers the paper publishes.
- [x] `paths.local.json` is untracked (`.gitignore`). No active machine-specific path defaults remain in code: path settings come from `paths.local.json` or environment variables. Explanatory references to the original `C:/SAT` locations remain in docstrings and provenance notes, and `paths.example.json` contains placeholder Windows paths. The Allen Coral Atlas path is now a setting (`aca_bathymetry`).
- [x] **Commit identity:** every commit after `main` uses the author's academic address. The three `main` commits are already on GitHub and are not rewritten.
- [x] **README:**
  - the SUBMITTED paper's results came from the author's original MATLAB workflow (surviving scripts and archived outputs listed in `docs/PROVENANCE.md`, not in this repository); the R1 paper's tables were rebuilt in Python by `r1/`; the earlier two-class MATLAB prototype was removed (history: commit 2c0d556, tag `matlab-prototype`);
  - Python is the verified reimplementation (`compat` reproduces the submitted Tables 4–5; `fixed` corrects the input scaling and is used by `r1/`);
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
- [x] `r1/` holds the R1 runner, `legacy_replay.py`, `table3_retest.py`, `table3_spatial_cv.py`, `make_tables.py`, `calibration_sensitivity.py`, `scenes.csv`/`scenes.py`, three configurations and `r1/README.md`. Configuration paths are placeholders; every SHA-256 and size is unchanged.
- [x] README results framing: R1 tables come from `r1/`; `compat` reproduces the submitted tables (provenance); `fixed` corrects the input scaling.
- [x] `docs/INPUTS.md` lists StableSand and the coregistered UDM2 flag files by hash and size.
- [x] Codex pre-push review of the `r1-release` diff (21 Nov 2026).
- [x] Every R1 table and supplement item is mapped to a release script in `r1/README.md` and was reproduced from the stored outputs (21 Nov 2026), except the archived S2 binary-analysis table.
- [x] Code identity hashes are line-ending neutral (LF).
- [x] Codex pre-push review (21 Nov 2026) applied: documented run order for Table 3 (`make_tables.py` into `<RUN>/report` before `table3_retest.py`); `tzdata` added for Windows; all `r1/` outputs pass `paths.assert_writable` and refuse configured input locations, including the failure log; exact-commit tagging with the CHANGELOG dated first; preflight wording (no scientific output before checks); submitted vs R1 provenance in this checklist; `validate` command for staged runs.

## Push (author)
1. Push only after the R1 manuscript is final and the author gives the go-ahead.
2. **Connect the repository to Zenodo first** (zenodo.org, GitHub integration, switch the repository on). A release made before this mints no DOI.
3. **Date the release in the commit that will be tagged.** Replace "unreleased; date set when tagged" in `CHANGELOG.md` with the release date and commit it on `r1-release`. Record that commit: `REL=$(git rev-parse r1-release)`.
4. **Check the remote, then fast-forward `main`** (no force, no PR merge, so `main` becomes exactly `$REL`):

   ```bash
   git fetch github
   git merge-base --is-ancestor github/main "$REL" && echo "fast-forward OK"   # must print OK
   git push github "$REL:refs/heads/main"
   git fetch github && test "$(git rev-parse github/main)" = "$REL" && echo "main = $REL"
   ```

5. **Tag that exact commit** and push the tags:

   ```bash
   git tag -a v1.0.0 "$REL" -m "CoralRemote 1.0.0 (AJSTR R1)"
   git push github v1.0.0
   git push github matlab-prototype   # annotated tag on 2c0d556 (created locally): link to the removed MATLAB prototype
   git ls-remote --tags github | grep -E 'v1\.0\.0\^\{\}|matlab-prototype\^\{\}'   # peeled tags: must show $REL and 2c0d556
   ```

   Then create the GitHub release from `v1.0.0`; Zenodo archives that tag and mints the DOI.
6. Add the DOI badge to the README (one more commit and push; the archived v1.0.0 stays as tagged). Add the repository URL and DOI to the paper's Data Availability statement.
