# Evaluation Scope

This folder is intentionally minimal in the first public repository.

The original private workspace contained A3-specific evaluation and experiment helpers here, but they were excluded from the public candidate because they were tightly coupled to private machine paths, private run artifacts, and local data products.

For the first GitHub release, the production methodology path is:

- `01_preprocessing`
- `02_baseline_construction`
- `03_classification`

Evaluation code can be reintroduced later after it is refactored to use repository config rather than private local paths.
