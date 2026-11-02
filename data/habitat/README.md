# Appendix C habitat inputs

These files are the compact, publication-ready spatial inputs used to
reproduce Appendix C of the manuscript.

- `camera_locations.gpkg`: 60 camera stations, with camera, sampling-polygon,
  and management-category identifiers. Two blank `POINT EMPTY` records and
  machine-specific path fields in the source shapefile were excluded.
- `sampling_polygons.gpkg`: the 12 sampling polygons, transformed to EPSG:2039
  to match the land-cover raster.
- `maarag_landcover_25m_study_extent.tif`: study-extent crop of the finest
  available categorical Maarag land-cover raster. The original 25 m cells,
  raster origin, class codes, and EPSG:2039 coordinate reference system are
  unchanged. The raster was cropped only; it was not reprojected or resampled.
- `landcover_key.csv`: English class labels and the explicit structural-cover
  definitions used in Appendix C. Structural woody cover is codes 33, 34, and
  35. Natural open/herbaceous cover is codes 31 and 32. Code -128 is NoData.
- `field_vegetation_structure.csv`: one field-assessed vegetation-structure
  classification per camera, derived from `ObsCovs_050524.csv` after verifying
  that each camera had one consistent non-missing classification.
- `source_provenance.csv`: checksums connecting these compact inputs to the
  original GIS and covariate files used to prepare them.

The crop extends beyond every 1,000 m camera buffer used in the sensitivity
analysis. Exact raster-cell intersection fractions are calculated by
`R/habitat_appendix.R`; class proportions use valid raster area as the
denominator.

Copyright (c) 2026 the manuscript authors. All rights reserved; see the
repository `LICENSE`.
