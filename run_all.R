required_packages <- c(
  "activity",
  "car",
  "dplyr",
  "emmeans",
  "ggplot2",
  "ggtext",
  "glmmTMB",
  "janitor",
  "overlap",
  "patchwork",
  "purrr",
  "readr",
  "Rcpp",
  "scales",
  "tibble",
  "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following R packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

source(file.path("R", "00_config.R"))
source(file.path("R", "01_helpers.R"))
source(file.path("R", "02_data.R"))
source(file.path("R", "part1_land_use_activity.R"))
source(file.path("R", "part2_nocturnality.R"))
source(file.path("R", "part4_human_activity.R"))
source(file.path("R", "part5_species_overlap.R"))
source(file.path("R", "part6_gazelle_focus.R"))
source(file.path("R", "build_manuscript_outputs.R"))
source(file.path("R", "validate_outputs.R"))

set.seed(SEED)
message("Golan diel-activity manuscript pipeline")
message("Run profile: ", RUN_PROFILE)
message("Output directory: ", OUTPUT_ROOT)

analysis_data <- load_analysis_data()
part1_results <- run_part1(analysis_data)
part2_results <- run_part2(analysis_data)
part4_results <- run_part4(analysis_data)
part5_results <- run_part5(analysis_data)
part6_results <- run_part6(analysis_data, part4_results)
detection_results <- build_detection_outputs(analysis_data)
table_b1_results <- build_table_b1(analysis_data)

results <- list(
  data = analysis_data,
  part1 = part1_results,
  part2 = part2_results,
  part4 = part4_results,
  part5 = part5_results,
  part6 = part6_results,
  detections = detection_results,
  table_b1 = table_b1_results
)

validate_analysis(results)
write_reproducibility_metadata(results)

message("Completed successfully: ", normalizePath(OUTPUT_ROOT))
