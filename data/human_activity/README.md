# Human-activity index inputs

These files are the archived inputs used to derive the camera-level
`human.abundance` covariate:

- `Human_Abundance_Normalized.csv`: repeated independent human-detection
  counts for 60 cameras and seven sampling occasions.
- `ObsCovs_050524.csv`: occasion-level detection covariates in camera-by-round
  order.
- `SiteCovs_060324.csv`: site covariates and camera identifiers used when the
  N-mixture model was fitted.

`R/human_activity_index.R` reconstructs the final Poisson N-mixture model from
the archived fitted coefficients, extracts empirical-Bayes site estimates,
scales them by their maximum, and verifies the result against `human.abundance`
in `data/SiteCovs_280825.csv`. A full likelihood refit is available through
`Rscript run_human_activity_index.R --refit-model`. The final fitted detection
formula is:

```r
~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1
```

The source files are retained without modification. Their MD5 checksums are
validated by the reconstruction script. The archived candidate-model table is
retained as
`manuscript_outputs/additional_results/human_detection_candidate_models_archived.csv`;
its MD5 is `c1e5597c8fd2568f87428cb6ffca78a7`.
