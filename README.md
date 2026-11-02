# Mammalian temporal niches and activity shift under culling and land management

Reproducible code, input data, manuscript figures, and statistical tables for
the submitted manuscript:

Shlomo Preiss-Bloom, Hila Shamon, Zachary Amir, Eamonn I.F. Wooster,
Dror Ben-Ami, and Tamar Dayan.

This repository accompanies a manuscript currently under review. It does not
imply acceptance or publication by any journal.

## Repository contents

- `data/`: canonical camera-event, site-covariate, Appendix B human-activity,
  and Appendix C habitat inputs.
- `R/`: modular analysis and output-generation functions.
- `src/`: compiled backend used to accelerate kernel randomization tests.
- `run_all.R`: complete analysis entry point.
- `manuscript_outputs/`: publication-profile figures, tables, complete results,
  diagnostics, and reproducibility metadata.

Part 3 peak analyses and post-submission exploratory analyses are not included
because they are not part of the submitted study. The scripts in `R/` are the
canonical analysis; exploratory and legacy scripts are intentionally excluded
from the publication repository.

## Reproduce the manuscript outputs

Clone the repository, open a terminal in its root, and restore the recorded R
environment:

```r
install.packages("renv")
renv::restore()
```

Then run:

```powershell
Rscript run_all.R
```

For a clean-checkout reproducibility audit, confirm that all analytical checks
pass and that every non-PDF output matches the release exactly:

```powershell
Rscript -e "x <- read.csv('manuscript_outputs/diagnostics/validation_checks.csv'); stopifnot(nrow(x) > 0, all(x[, 'passed']))"
git diff --exit-code -- manuscript_outputs ":(exclude)**/*.pdf" ":(exclude)manuscript_outputs/metadata/output_manifest.csv"
```

The Cairo PDF device may write different internal compressed metadata between
runs, and `output_manifest.csv` consequently may report slightly different PDF
byte sizes. These binary differences do not affect the plotted content. The
corresponding PNG figures, statistical CSVs, validation records, and other
non-PDF outputs must remain unchanged.

The default `publication` profile uses the manuscript settings: 9,999 Part 1
randomizations, 999 Part 4 randomizations, 9,999 Part 6 randomizations, 4,999
nocturnality bootstraps, 499 interspecific-overlap bootstraps, and 1,999
Delta-hat 1 uncertainty bootstraps. Low-event Delta-hat 1 tests use 999
permutations. Set `GOLAN_CORES` to a positive integer to control parallel
workers; the default is at most four workers to avoid repeated Windows PSOCK
connection exhaustion.

For a fast installation and code-path check only:

```powershell
$env:GOLAN_RUN_PROFILE = "check"
Rscript run_all.R
Remove-Item Env:GOLAN_RUN_PROFILE
```

Check-profile output is not suitable for inference or publication.

Appendix B's site-level human-activity index can be reconstructed separately:

```powershell
Rscript run_human_activity_index.R
```

This reconstructs the index from the archived fitted coefficients and verifies
the scaled site estimates against `human.abundance` in `SiteCovs_280825.csv`.
The selected likelihood can be fully refitted with the same model and data by
running `Rscript run_human_activity_index.R --refit-model`. Refitting every
candidate detection model is also optional because those likelihood fits are
slow:

```powershell
Rscript run_human_activity_index.R --refit-candidates
```

## Analytical conventions

Camera timestamps are camera-local `Asia/Jerusalem` times and are not converted
from UTC. The pipeline checks the canonical input checksum and a known Camelot
timestamp to prevent recurrence of the previously identified time-zone shift.

Exact camera-species-timestamp records repeated across adjacent Camelot session
IDs are collapsed before focal-species analyses. This applies the manuscript's
30-minute independent-event definition consistently across all analyses. The
collapsed source rows are recorded in
`manuscript_outputs/diagnostics/focal_exact_cross_session_duplicates_collapsed.csv`.

Figure 2 and descriptive Table A.2 retain the source-export total of 4,334 wild
boar detections. Inferential focal-species analyses collapse the duplicated
cross-session row and therefore use 4,333 independent wild boar events. This
distinction is asserted by the validation pipeline and documented in the
diagnostic table.

Activity curves are circular probability-density functions normalized to area
one and weighted by observed group size. Delta-hat 1 is used when the smaller
event sample has fewer than 50 observations; Delta-hat 4 is used otherwise.
Nocturnality is absolute circular distance in hours from solar noon, where 0 is
solar noon and 12 is solar midnight.

For descriptive summaries, land-use means and 95% bootstrap CIs are estimated
by resampling independent detections with replacement and recalculating
group-size-weighted means. Figure A.1 applies the same procedure independently
to each land-use category and its minefield reference, then calculates each
replicate contrast. The metadata and table captions distinguish the uncertainty
procedures used by other analyses.

Appendix C uses the actual camera coordinates and exact fractional extraction
from the 25 m categorical Maarag raster at 250, 500, and 1,000 m. The 500 m
scale is primary. Camera values are averaged within each sampling polygon, and
the 12 polygons are the independent units in exact tests that enumerate all
369,600 assignments to four labelled groups of three. Structural woody cover
comprises shrubs, Mediterranean maquis, and planted conifer forest; orchards
are excluded. The public raster is a study-extent crop of the source raster,
with no reprojection or resampling.

Appendix B's human-activity index uses a Poisson N-mixture model with an
intercept-only state process and a detection model containing path type,
camera model, vegetation structure, and passability behind the camera. The
recovered original count history and covariates are checksum-protected, and
the reconstruction script validates the resulting scaled index against the
covariate used by all downstream analyses.

## Generated outputs

`manuscript_outputs/` contains:

- `main_figures/`: Figures 1-6.
- `main_tables/`: Table 1.
- `supplementary_figures/`: Figures A.1-A.2.
- `supplementary_tables/`: Tables A.1-A.5, B.1, and Appendix C habitat
  summaries.
- `additional_results/`: complete, unrounded statistical results.
- `diagnostics/`: sample assignments, density-area audits, and validation.
- `metadata/`: settings, checksums, package versions, and output manifest.

Figure 1 is an external QGIS cartographic product. The included JPEG is the
version used with the manuscript; its provenance and checksum are recorded in
`manuscript_outputs/metadata/external_assets.csv`.

## Citation

Until the manuscript receives a permanent bibliographic record, cite this
repository using `CITATION.cff` and the archived software release. The citation
metadata can be updated with the article DOI after acceptance.

## License

Copyright (c) 2026 the manuscript authors. All rights reserved. See `LICENSE`.
