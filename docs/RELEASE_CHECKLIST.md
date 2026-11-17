# Release checklist for v1.0.0

**Status (17 Nov 2026): ready to push; nothing has been pushed, tagged or published.**
- The release branch is `python-port`. It fast-forwards from GitHub `main` (2c0d556), so no history is rewritten on GitHub.
- The `github` remote is configured for fetching. Pushing needs the author's go-ahead.

## Done
- [x] **No data tracked.** `git ls-files | grep -Ei '\.(mat|tif|tiff|npy|npz|h5)$'` prints nothing. `tests/golden/regression_tile.json` holds only aggregates (400-pixel means, fractions, class counts), the same kind of numbers the paper publishes.
- [x] `paths.local.json` is untracked (`.gitignore`). No active machine-specific path defaults remain in code: path settings come from `paths.local.json` or environment variables. Explanatory references to the original `C:/SAT` locations remain in docstrings and provenance notes, and `paths.example.json` contains placeholder Windows paths. The Allen Coral Atlas path is now a setting (`aca_bathymetry`).
- [x] **Commit identity:** every commit after `main` uses the author's academic address. The three `main` commits are already on GitHub and are not rewritten.
- [x] **README:**
  - the paper's results came from the author's original MATLAB workflow (surviving scripts and archived outputs listed in `docs/PROVENANCE.md`, not in this repository); `legacy/matlab/` is an earlier public two-class workflow, not the paper's code;
  - Python is the verified reimplementation (`compat` reproduces Tables 4–5; `fixed` corrects the input scaling);
  - the units numbers match the paper's paired comparison (§3.5.6).
- [x] **Version 1.0.0:** set in `pyproject.toml` and `coralremote/__init__.py`. The CHANGELOG is updated, and `CITATION.cff` has been added.
- [x] **Internal port reports kept** (`docs/PY4–PY8`, `R1_LOO_DATE`). Local paths were removed. Their R1 details are the ones the revised paper discloses.
- [x] **No second-manuscript material:** the change-detection work lives in a separate local repository, `coralremote-change`, which is not part of this release.

## Push (author)
1. Push only after the R1 manuscript is final, as decided on 17 Nov.
2. Push the release branch:

   ```bash
   git push github python-port
   ```

   Then open a PR into `main` on GitHub and merge it. Alternatively, fast-forward `main` directly: `git push github python-port:main`.
3. Before creating a release, connect the repository to Zenodo (GitHub integration) so a DOI is minted.
4. Tag and release:

   ```bash
   git tag -a v1.0.0 -m "CoralRemote 1.0.0"
   git push github v1.0.0
   ```

   Then create the GitHub release from the tag and date the CHANGELOG entry.
5. Add the DOI badge to the README. At acceptance (or in round 2), replace "available from the corresponding author upon reasonable request" in the paper with the repository URL and DOI.
