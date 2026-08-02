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

The source export contains one focal-species event repeated at the same camera
and timestamp across adjacent Camelot session IDs. The shared data loader
collapses exact cross-session repeats consistently before all focal-species
analyses and records the affected source rows in the diagnostics output.

## SiteCovs_280825.csv

Camera-level land-use assignments and standardized activity or abundance
covariates used to define pressure gradients and Q1/Q4 camera groups.

The dataset contains 60 unique camera stations, with 15 cameras in each of four
land-use categories.

## Data availability and rights

The files include exact camera coordinates and event-level observations. Their
public release has been authorized by the authors. Copyright (c) 2026 the
manuscript authors. All rights reserved; see the repository `LICENSE`.
