# Release checklist for v1.0.0 (author actions)

Everything below needs the author's explicit OK. Nothing has been pushed, tagged or published. Branch `python-port` lives only in the local clone, and its push URL is disabled.

## Before pushing
- [ ] Review the diff of `python-port` against `main`:
  - the MATLAB tree moved to `legacy/matlab/`;
  - the new `coralremote/`, `tests/`, `scripts/` and `docs/`.
- [ ] Decide whether the internal port reports (`docs/PY4_*.md`, `PY5_*`, `PY6_*`) are public. They mention local folder names and unpublished R1 details. Options: keep them, trim them, or move them out of the repo.
- [ ] Confirm that no data is tracked:

  ```bash
  git ls-files | grep -Ei '\.(mat|tif|tiff|npy|npz|h5)$'
  ```

  This should print nothing. `tests/golden/regression_tile.json` holds only aggregates (means, fractions, class counts). Confirm that's acceptable under the Planet licence.
- [ ] Check that `paths.local.json` is untracked (it's in `.gitignore`).
- [ ] Decide on the AJSTR paper (plan §5): the units fix is material under the >2 pp Table 5 rule. Should the repository state this publicly at release, or only after the editor has been told?

## Publish
1. Re-enable the push URL and push `python-port`; open a PR into `main` and merge after review.
2. Set `version = "1.0.0"` in `pyproject.toml` and `coralremote/__init__.py`, and date the CHANGELOG entry.
3. Tag and create the GitHub release:

   ```bash
   git tag -a v1.0.0 -m "CoralRemote 1.0.0"
   ```

4. Connect the repository to Zenodo (GitHub integration) **before** creating the release, so a DOI is minted. Add the DOI badge and citation to the README.
5. Optionally add a `CITATION.cff`.
