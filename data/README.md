# Input data

This directory contains the two canonical inputs used by `run_all.R`.

## Independent_Captures_Sep_2025_solar.csv

Camera-trap event records with camera ID, Camelot deployment-session ID,
species identity, camera-local timestamp, observed group size, coordinates,
and solar-standardized circular time.

The timestamps were recorded as local `Asia/Jerusalem` camera time. They must
not be interpreted as UTC or shifted by an additional UTC offset. The analysis
checks the file MD5 (`f54e92df14a6bcb1730ecb5474370db4`) and a known timestamp
before proceeding.

The source export contains one wild boar event repeated at the same camera and
timestamp across adjacent Camelot session IDs. The shared data loader collapses
that repeat before inferential focal-species analyses and records both source
rows in the diagnostics output. Figure 2 and descriptive Table A.2 retain the
source-export total (4,334); inferential analyses use 4,333 independent events.

## SiteCovs_280825.csv

Camera-level land-use assignments and standardized activity or abundance
covariates used to define pressure gradients and Q1/Q4 camera groups.

The dataset contains 60 unique camera stations, with 15 cameras in each of four
land-use categories.

The `human.abundance` column is the N-mixture covariate used to rank cameras.
The repeated-count history, observation covariates, and reconstruction code
are provided in `human_activity/` and `R/human_activity_index.R`. The script
verifies the reconstructed index against this column before returning.

## human_activity/

Inputs used to derive the site-level human-activity index described in
Appendix B. Run `Rscript run_human_activity_index.R` after restoring the `renv`
environment. Add `--refit-model` to rerun the final likelihood optimization.
Add `--refit-candidates` to refit the full candidate model set; both options
are computationally intensive.

## habitat/

Publication-ready spatial inputs for Appendix C. See
`data/habitat/README.md` for class definitions, provenance, and the documented
study-extent crop of the original 25 m Maarag raster.

## Data availability and rights

The files include exact camera coordinates and event-level observations. Their
public release has been authorized by the authors. Copyright (c) 2026 the
manuscript authors. All rights reserved; see the repository `LICENSE`.
