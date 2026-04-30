# Code Section Map

This `code` folder is organized into four cautious publish sections.

## Sections

- `01_preprocessing`
- `02_baseline_construction`
- `03_classification`
- `04_evaluation`
- `utils`

## Main entry points

- `RunAll_LogRatio.m`
- `RunAll_LogRatio_A3.m`
- `RunAll_LogRatio_Mode.m`
- `01_preprocessing/Preprocess_Run.m`
- `02_baseline_construction/Build_Baseline.m`
- `02_baseline_construction/Build_BaselineCoralMask.m`
- `03_classification/Train_Stage1.m`
- `03_classification/Classify_Stage1.m`
- `03_classification/Train_Stage2.m`
- `03_classification/Classify_Stage2.m`
- `03_classification/Analyze_Change.m`

## Preprocessing

The public preprocessing entry point is `01_preprocessing/Preprocess_Run.m`.
It keeps the tested step scripts intact and runs them in this order:

- `Preprocess_Masking.m`
- `Preprocess_Normalize.m`
- `Preprocess_WaveMask.m`
- `Preprocess_FeatureStack.m`

## Scope

This public candidate keeps the methodology pipeline first.
The runnable public path is concentrated in:

- `01_preprocessing`
- `02_baseline_construction`
- `03_classification`

`04_evaluation` is intentionally trimmed in the first release and now contains only a scope note.
