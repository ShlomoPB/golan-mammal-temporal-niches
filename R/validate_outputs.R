validate_analysis <- function(results) {
  message("Validating analytical and graphical outputs")
  checks <- tribble(
    ~check, ~passed, ~detail,
    "Seven focal species",
    n_distinct(results$data$focal$species) == 7,
    paste("observed", n_distinct(results$data$focal$species)),
    "Four land-use categories",
    n_distinct(results$data$focal$land_use) == 4,
    paste("observed", n_distinct(results$data$focal$land_use)),
    "Figure labels use the seven final common names",
    identical(
      unname(SPECIES_COMMON[SPECIES]),
      c(
        "Grey wolf",
        "Golden jackal",
        "Wild boar",
        "Mountain gazelle",
        "Red fox",
        "Caucasian badger",
        "Indian crested porcupine"
      )
    ),
    "common-name taxonomy and species order must match the revised manuscript",
    "60 camera stations",
    n_distinct(results$data$site_covs$camera_name) == 60,
    paste("observed", n_distinct(results$data$site_covs$camera_name)),
    "One exact cross-session focal event is collapsed",
    nrow(results$data$focal_exact_duplicate_rows) == 2 &&
      n_distinct(results$data$focal_exact_duplicate_rows$species) == 1 &&
      as.character(results$data$focal_exact_duplicate_rows$species[1]) ==
        "Sus scrofa",
    paste(
      nrow(results$data$focal_exact_duplicate_rows),
      "export rows representing one wild-boar event"
    ),
    "Table 1 has 42 species-by-pair rows",
    nrow(results$part1$table1) == 42,
    paste("observed", nrow(results$part1$table1)),
    "Table A.1 has seven species",
    nrow(results$part4$table_a1) == 7,
    paste("observed", nrow(results$part4$table_a1)),
    "Table A.2 has seven species plus humans",
    nrow(results$detections$table_a2) == 8,
    paste("observed", nrow(results$detections$table_a2)),
    "Table A.2 reports the submitted descriptive wild-boar total",
    with(
      results$detections$table_a2 |>
        filter(species == "Sus scrofa"),
      independent_observations == 4334
    ),
    "source export n = 4,334; inferential focal data n = 4,333",
    "Table A.3 has 28 species-by-land-use rows",
    nrow(results$part2$species_land_use) == 28,
    paste("observed", nrow(results$part2$species_land_use)),
    "Figure A.1 has 28 species-by-land-use contrast rows",
    nrow(results$part2$figure_a1_contrasts) == 28,
    paste("observed", nrow(results$part2$figure_a1_contrasts)),
    "Figure A.1 minefield references are zero without uncertainty bars",
    with(
      results$part2$figure_a1_contrasts |>
        filter(as.character(land_use) == "Minefields"),
      all(contrast_vs_minefields == 0) &&
        all(is.na(ci_low)) &&
        all(is.na(ci_high))
    ),
    "all seven minefield references must be fixed at zero",
    "Figure A.1 point estimates match corrected land-use means",
    all(abs(
      results$part2$figure_a1_contrasts$mean_nocturnality -
        results$part2$species_land_use$mean
    ) < 1e-12),
    "absolute tolerance = 1e-12",
    "Figure A.1 sample sizes match Table A.3",
    all(
      results$part2$figure_a1_contrasts$n_detections ==
        results$part2$species_land_use$n_events &
        results$part2$figure_a1_contrasts$sum_group_size_weights ==
          results$part2$species_land_use$n_individuals
    ),
    "event counts and sums of group-size weights must match",
    "Every Figure A.1 estimate lies within its contrast CI",
    RUN_PROFILE != "publication" || with(
      results$part2$figure_a1_contrasts |>
        filter(as.character(land_use) != "Minefields"),
      all(contrast_vs_minefields >= ci_low & contrast_vs_minefields <= ci_high)
    ),
    if (RUN_PROFILE == "publication") {
      "all 21 non-reference estimates must be bounded by their CIs"
    } else {
      "publication-only CI-bound check skipped for nine-draw check profile"
    },
    "Table A.4 source has 84 pair-by-land-use rows",
    nrow(results$part5$overlaps) == 84,
    paste("observed", nrow(results$part5$overlaps)),
    "All Part 1 density curves integrate to one",
    all(abs(
      results$part1$curves |>
        group_by(species, land_use) |>
        summarise(area = trapz(hour, density), .groups = "drop") |>
        pull(area) - 1
    ) < 0.002),
    "tolerance = 0.002",
    "All Part 4 density curves integrate to one",
    all(abs(
      results$part4$curves |>
        group_by(species, group) |>
        summarise(area = trapz(hour, density), .groups = "drop") |>
        pull(area) - 1
    ) < 0.002),
    "tolerance = 0.002",
    "No Table A.3 mean lies outside its 95% CI",
    RUN_PROFILE != "publication" || all(with(
      results$part2$species_land_use,
      mean >= ci_low & mean <= ci_high
    )),
    if (RUN_PROFILE == "publication") {
      "mean must be bounded by CI"
    } else {
      "publication-only inferential check skipped for nine-draw check profile"
    },
    "Gazelle Q1/Q4 groups contain exactly 15 cameras each",
    length(results$part4$human_groups$low) == 15 &&
      length(results$part4$human_groups$high) == 15,
    paste(
      length(results$part4$human_groups$low),
      "and",
      length(results$part4$human_groups$high)
    ),
    "Gazelle Part 4 event samples match the exact Q1/Q4 split",
    with(
      results$part4$human_summary |>
        filter(species == "Gazella gazella"),
      n_low_Q1 == 960 && n_high_Q4 == 483
    ),
    "expected Q1 = 960 and Q4 = 483 events",
    "Badger Part 4 event samples match the exact Q1/Q4 split",
    with(
      results$part4$human_summary |>
        filter(species == "Meles meles"),
      n_low_Q1 == 58 && n_high_Q4 == 49
    ),
    "expected Q1 = 58 and Q4 = 49 events",
    "Wild-boar Part 4 samples exclude the duplicated session-boundary event",
    with(
      results$part4$human_summary |>
        filter(species == "Sus scrofa"),
      n_low_Q1 == 932 && n_high_Q4 == 886
    ),
    "expected Q1 = 932 and Q4 = 886 independent events",
    "Gazelle Q4-Q1 nocturnality effect matches camera-local timestamps",
    abs(
      results$part6$gazelle_summary$delta_noct_Q4_minus_Q1 -
        1.067772246
    ) < 1e-6,
    sprintf(
      "Q4-Q1 = %.4f h",
      results$part6$gazelle_summary$delta_noct_Q4_minus_Q1
    ),
    "Gazelle overlap agrees with the shaded unit-area curves",
    abs(
      results$part6$gazelle_summary$overlap_delta -
        curve_overlap(
          results$part4$curves |>
            filter(
              species == "Gazella gazella",
              group == "Low human activity (Q1)"
            ) |>
            select(hour, density),
          results$part4$curves |>
            filter(
              species == "Gazella gazella",
              group == "High human activity (Q4)"
            ) |>
            select(hour, density)
        )
    ) < 0.01,
    "absolute tolerance = 0.01",
    "Table B.1 has 15 cameras in every quartile",
    all(results$table_b1$counts$Total == 15),
    paste(results$table_b1$counts$Total, collapse = ", "),
    "Appendix C habitat design has 60 cameras and 12 polygons",
    results$habitat$n_cameras == 60L && results$habitat$n_polygons == 12L,
    paste(
      results$habitat$n_cameras,
      "cameras and",
      results$habitat$n_polygons,
      "polygons"
    ),
    "Appendix C structural-woody test matches the submitted result",
    with(
      results$habitat$tests |>
        filter(buffer_radius_m == 500L, variable == "structural_woody"),
      abs(f_statistic - 6.999404287425887) < 1e-10 &&
        abs(eta_squared - 0.7241209298657573) < 1e-10 &&
        abs(p_exact - 0.026623376623376622) < 1e-12
    ),
    "expected F = 6.999, eta-squared = 0.724, P = 0.027",
    "Appendix C natural-open test matches the submitted result",
    with(
      results$habitat$tests |>
        filter(buffer_radius_m == 500L, variable == "natural_open"),
      abs(f_statistic - 4.860524921507623) < 1e-10 &&
        abs(eta_squared - 0.6457288703988651) < 1e-10 &&
        abs(p_exact - 0.04896103896103896) < 1e-12
    ),
    "expected F = 4.861, eta-squared = 0.646, P = 0.049",
    "Appendix C field-vegetation test matches the submitted result",
    abs(results$habitat$field_test$p_exact - 0.5761038961038961) < 1e-12,
    "expected exact P = 0.576",
    "All Part 5 beta models report optimizer convergence",
    all(results$part5$model_diagnostics$optimizer_convergence_code == 0),
    paste(
      results$part5$model_diagnostics$optimizer_convergence_code,
      collapse = ", "
    ),
    "All Part 5 beta models have positive-definite Hessians",
    all(results$part5$model_diagnostics$positive_definite_hessian),
    paste(
      results$part5$model_diagnostics$positive_definite_hessian,
      collapse = ", "
    )
  )

  expected_files <- c(
    file.path(OUTPUT_DIRS$main_figures, "Figure_2_detection_totals_by_land_use.png"),
    file.path(OUTPUT_DIRS$main_figures, "Figure_3_land_use_activity.png"),
    file.path(OUTPUT_DIRS$main_figures, "Figure_4_nocturnality_gradient.png"),
    file.path(OUTPUT_DIRS$main_figures, "Figure_5_pairwise_overlap_matrices.png"),
    file.path(
      OUTPUT_DIRS$main_figures,
      "Figure_6_gazelle_human_and_reserve_activity.png"
    ),
    file.path(OUTPUT_DIRS$main_tables, "Table_1_land_use_activity_tests.csv"),
    file.path(
      OUTPUT_DIRS$supplementary_figures,
      "Figure_A1_nocturnality_by_land_use.png"
    ),
    file.path(
      OUTPUT_DIRS$supplementary_figures,
      "Figure_A2_activity_by_human_activity.png"
    ),
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_C1_habitat_composition_500m.csv"
    )
  )
  checks <- bind_rows(
    checks,
    tibble(
      check = paste("Output exists:", basename(expected_files)),
      passed = file.exists(expected_files),
      detail = {
        portable_paths <- normalizePath(
          expected_files,
          winslash = "/",
          mustWork = FALSE
        )
        project_prefix <- paste0(
          normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE),
          "/"
        )
        ifelse(
          startsWith(portable_paths, project_prefix),
          substring(portable_paths, nchar(project_prefix) + 1L),
          portable_paths
        )
      }
    )
  )

  write_csv(checks, file.path(OUTPUT_DIRS$diagnostics, "validation_checks.csv"))
  failed <- checks |> filter(!passed)
  if (nrow(failed)) {
    print(failed, n = Inf)
    stop("Validation failed; see diagnostics/validation_checks.csv.")
  }

  invisible(checks)
}
