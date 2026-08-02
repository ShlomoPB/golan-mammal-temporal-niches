build_detection_outputs <- function(data) {
  message("Building Figure 2 and Table A.2")
  captures <- data$captures
  focal <- data$focal
  site_covs <- data$site_covs

  focal_summary <- map_dfr(SPECIES, function(sp) {
    species_data <- focal |> filter(species == sp)
    total <- sum(species_data$sighting_quantity)
    by_land_use <- species_data |>
      group_by(land_use) |>
      summarise(n = sum(sighting_quantity), .groups = "drop") |>
      complete(
        land_use = factor(LAND_USE_LEVELS, levels = LAND_USE_LEVELS),
        fill = list(n = 0)
      )
    tibble(
      species = sp,
      common_name = unname(SPECIES_COMMON[sp]),
      proportion_sites = n_distinct(species_data$camera_name) / nrow(site_covs),
      independent_observations = total,
      High_culling = by_land_use$n[by_land_use$land_use == "High culling"] / total,
      Low_culling = by_land_use$n[by_land_use$land_use == "Low culling"] / total,
      Nature_reserves = by_land_use$n[by_land_use$land_use == "Nature reserves"] / total,
      Minefields = by_land_use$n[by_land_use$land_use == "Minefields"] / total
    )
  })

  humans <- captures |>
    filter(
      species %in% HUMAN_ALIASES,
      !is.na(land_use),
      is.finite(sighting_quantity),
      sighting_quantity > 0
    )
  human_total <- sum(humans$sighting_quantity)
  human_by_land_use <- humans |>
    group_by(land_use) |>
    summarise(n = sum(sighting_quantity), .groups = "drop") |>
    complete(
      land_use = factor(LAND_USE_LEVELS, levels = LAND_USE_LEVELS),
      fill = list(n = 0)
    )
  human_summary <- tibble(
    species = "Homo sapiens",
    common_name = "Human",
    proportion_sites = n_distinct(humans$camera_name) / nrow(site_covs),
    independent_observations = human_total,
    High_culling =
      human_by_land_use$n[human_by_land_use$land_use == "High culling"] / human_total,
    Low_culling =
      human_by_land_use$n[human_by_land_use$land_use == "Low culling"] / human_total,
    Nature_reserves =
      human_by_land_use$n[human_by_land_use$land_use == "Nature reserves"] / human_total,
    Minefields =
      human_by_land_use$n[human_by_land_use$land_use == "Minefields"] / human_total
  )

  table_a2 <- bind_rows(focal_summary, human_summary)
  write_csv(
    table_a2,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_A2_detection_summary_and_naive_occupancy.csv"
    )
  )

  figure_data <- focal_summary |>
    pivot_longer(
      cols = c(High_culling, Low_culling, Nature_reserves, Minefields),
      names_to = "land_use",
      values_to = "proportion"
    ) |>
    mutate(
      land_use = recode(
        land_use,
        High_culling = "High culling",
        Low_culling = "Low culling",
        Nature_reserves = "Nature reserves",
        Minefields = "Minefields"
      ),
      land_use = factor(land_use, levels = LAND_USE_LEVELS),
      detection_count = independent_observations * proportion,
      species = factor(species, levels = rev(DETECTION_SPECIES_ORDER)),
      common_name_plot = if_else(
        as.character(species) == "Hystrix indica",
        "Indian crested<br>porcupine",
        common_name
      ),
      label = paste0(
        "<b>", common_name_plot, "</b><br><i>", as.character(species), "</i>"
      )
    )
  label_levels <- focal_summary |>
    mutate(
      species = factor(species, levels = rev(DETECTION_SPECIES_ORDER)),
      common_name_plot = if_else(
        as.character(species) == "Hystrix indica",
        "Indian crested<br>porcupine",
        common_name
      ),
      label = paste0(
        "<b>", common_name_plot, "</b><br><i>", as.character(species), "</i>"
      )
    ) |>
    arrange(species) |>
    pull(label)
  figure_data$label <- factor(figure_data$label, levels = label_levels)

  label_data <- focal_summary |>
    mutate(
      common_name_plot = if_else(
        species == "Hystrix indica",
        "Indian crested<br>porcupine",
        common_name
      ),
      label = paste0("<b>", common_name_plot, "</b><br><i>", species, "</i>"),
      label = factor(label, levels = label_levels),
      right_label = paste0(
        "n = ", comma(independent_observations),
        "\nnaive occ. = ", sprintf("%.2f", proportion_sites)
      )
    )

  figure2 <- ggplot(
    figure_data,
    aes(detection_count, label, fill = land_use)
  ) +
    geom_col(
      width = 0.72,
      color = "white",
      linewidth = 0.35,
      position = position_stack(reverse = TRUE)
    ) +
    geom_text(
      data = label_data,
      aes(
        x = independent_observations + 350,
        y = label,
        label = right_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4.5,
      lineheight = 0.96
    ) +
    scale_fill_manual(
      values = DETECTION_FIGURE_COLORS,
      name = "Land-use category"
    ) +
    scale_x_continuous(
      breaks = seq(0, 10000, 2500),
      labels = comma,
      limits = c(0, 14200),
      expand = c(0, 0)
    ) +
    labs(x = "Independent detections", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 15),
      legend.text = element_text(size = 14),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = ggtext::element_markdown(
        color = "black",
        size = 13.5,
        lineheight = 0.95
      ),
      axis.title.x = element_text(size = 15),
      axis.text.x = element_text(size = 13),
      plot.margin = margin(8, 55, 8, 8)
    )
  save_figure(
    figure2,
    file.path(
      OUTPUT_DIRS$main_figures,
      "Figure_2_detection_totals_by_land_use"
    ),
    width = 10.8,
    height = 6.5
  )

  list(table_a2 = table_a2, figure2 = figure2)
}

build_table_b1 <- function(data) {
  message("Building Appendix B Table B.1")
  site_covs <- data$site_covs |>
    arrange(human.abundance, as.numeric(camera_name)) |>
    mutate(
      human_activity_quartile = factor(
        rep(c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"), each = 15),
        levels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")
      )
    )

  counts <- site_covs |>
    count(human_activity_quartile, land_use, name = "n_cameras") |>
    complete(
      human_activity_quartile,
      land_use = factor(LAND_USE_LEVELS, levels = LAND_USE_LEVELS),
      fill = list(n_cameras = 0)
    )
  count_wide <- counts |>
    pivot_wider(names_from = land_use, values_from = n_cameras) |>
    mutate(Total = rowSums(across(all_of(LAND_USE_LEVELS))))

  contingency <- xtabs(
    n_cameras ~ human_activity_quartile + land_use,
    data = counts
  )
  fisher <- fisher.test(contingency)
  test_table <- tibble(
    Test = "Fisher-Freeman-Halton exact test",
    p_value = fisher$p.value
  )

  write_csv(
    count_wide,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_B1_human_activity_quartile_by_land_use.csv"
    )
  )
  write_csv(
    test_table,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_B1_fisher_exact_test.csv"
    )
  )
  write_csv(
    site_covs |>
      select(
        camera_name,
        human.abundance,
        human_activity_quartile,
        land_use
      ),
    file.path(
      OUTPUT_DIRS$diagnostics,
      "Table_B1_camera_assignments.csv"
    )
  )

  list(counts = count_wide, test = test_table)
}

write_reproducibility_metadata <- function(results) {
  settings <- tibble(
    setting = c(
      "run_profile",
      "capture_time_basis",
      "capture_input_md5",
      "seed",
      "part1_compareCkern_randomizations",
      "part4_compareCkern_randomizations",
      "part6_compareCkern_randomizations",
      "Dhat1_permutations",
      "nocturnality_bootstrap_replicates",
      "interspecific_overlap_bootstrap_replicates",
      "Dhat1_overlap_bootstrap_replicates",
      "delta1_threshold_smaller_event_sample",
      "activity_curve_grid_points",
      "Q1_Q4_cameras_per_group",
      "parallel_workers"
    ),
    value = c(
      RUN_PROFILE,
      CAPTURE_TIME_BASIS,
      EXPECTED_CAPTURES_MD5,
      SEED,
      N_RANDOMIZATIONS_PART1,
      N_RANDOMIZATIONS_PART4,
      N_RANDOMIZATIONS_PART6,
      N_RANDOMIZATIONS_DHAT1,
      N_BOOT_NOCTURNALITY,
      N_BOOT_INTERSPECIES_OVERLAP,
      N_BOOT_DHAT1,
      DELTA1_THRESHOLD,
      N_ACTIVITY_GRID,
      N_Q_EXTREME_CAMERAS,
      N_CORES
    )
  )
  write_csv(settings, file.path(OUTPUT_DIRS$metadata, "run_settings.csv"))

  write_csv(
    results$data$focal_exact_duplicate_rows,
    file.path(
      OUTPUT_DIRS$diagnostics,
      "focal_exact_cross_session_duplicates_collapsed.csv"
    )
  )

  write_csv(
    tibble(
      backend = if (.fast_ckern_loaded) {
        "Rcpp/OpenMP von Mises kernel backend"
      } else {
        "activity package R fallback"
      },
      equivalence_tolerance = 1e-10,
      maximum_observed_equivalence_error =
        .fast_ckern_max_validation_error
    ),
    file.path(
      OUTPUT_DIRS$metadata,
      "compareCkern_backend_validation.csv"
    )
  )

  package_names <- c(
    "R",
    "activity",
    "overlap",
    "dplyr",
    "tidyr",
    "purrr",
    "readr",
    "janitor",
    "lubridate",
    "ggplot2",
    "ggtext",
    "patchwork",
    "glmmTMB",
    "emmeans",
    "car",
    "scales"
  )
  package_versions <- tibble(
    package = package_names,
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      map_chr(package_names[-1], ~as.character(packageVersion(.x)))
    )
  )
  write_csv(
    package_versions,
    file.path(OUTPUT_DIRS$metadata, "package_versions.csv")
  )

  external_figure <- file.path(
    OUTPUT_DIRS$main_figures,
    "Figure_1_map_external.jpeg"
  )
  if (file.exists(external_figure)) {
    write_csv(
      tibble(
        file = "main_figures/Figure_1_map_external.jpeg",
        provenance = paste(
          "Exact embedded Figure 1 extracted from the submitted manuscript;",
          "the QGIS project and source spatial layers are not in this repository."
        ),
        md5 = unname(tools::md5sum(external_figure))
      ),
      file.path(OUTPUT_DIRS$metadata, "external_assets.csv")
    )
  }

  environment_notes <- character()
  if (
    requireNamespace("glmmTMB", quietly = TRUE) &&
      requireNamespace("TMB", quietly = TRUE)
  ) {
    tmb_required <- tryCatch(
      as.character(
        get(
          ".TMB.build.version",
          envir = asNamespace("glmmTMB")
        )
      ),
      error = function(e) NA_character_
    )
    tmb_runtime <- as.character(utils::packageVersion("TMB"))
    environment_notes <- c(
      environment_notes,
      paste0(
        "glmmTMB runtime check: glmmTMB ",
        as.character(utils::packageVersion("glmmTMB")),
        "; TMB ",
        tmb_runtime,
        "; TMB version used to build glmmTMB: ",
        tmb_required,
        "; versions match: ",
        identical(tmb_required, tmb_runtime),
        "."
      )
    )
  }
  writeLines(
    environment_notes,
    file.path(OUTPUT_DIRS$metadata, "environment_notes.txt"),
    useBytes = TRUE
  )

  checksum_files <- c(
    DATA_CAPTURES,
    DATA_SITECOVS,
    file.path(PROJECT_DIR, "run_all.R"),
    file.path(PROJECT_DIR, "README.md"),
    file.path(PROJECT_DIR, "renv.lock"),
    file.path(PROJECT_DIR, "data", "README.md"),
    file.path(PROJECT_DIR, "DESCRIPTION"),
    file.path(PROJECT_DIR, "LICENSE"),
    file.path(PROJECT_DIR, "CITATION.cff"),
    list.files(
      file.path(PROJECT_DIR, "R"),
      pattern = "\\.R$",
      full.names = TRUE
    ),
    list.files(
      file.path(PROJECT_DIR, "src"),
      pattern = "\\.(cpp|h)$",
      full.names = TRUE
    )
  )
  checksum_files <- checksum_files[file.exists(checksum_files)]
  write_csv(
    tibble(
      file = sub(
        paste0(
          "^",
          gsub("([.()])", "\\\\\\1", PROJECT_DIR),
          "/?"
        ),
        "",
        normalizePath(checksum_files, winslash = "/", mustWork = TRUE)
      ),
      md5 = unname(tools::md5sum(checksum_files))
    ),
    file.path(OUTPUT_DIRS$metadata, "source_checksums.csv")
  )

  capture.output(
    sessionInfo(),
    file = file.path(OUTPUT_DIRS$metadata, "sessionInfo.txt")
  )

  generated <- list.files(
    OUTPUT_ROOT,
    recursive = TRUE,
    full.names = TRUE
  )
  manifest <- tibble(
    file = sub(
      paste0("^", gsub("([.()])", "\\\\\\1", OUTPUT_ROOT), "[/\\\\]?"),
      "",
      normalizePath(generated, winslash = "/", mustWork = TRUE)
    ),
    bytes = file.info(generated)$size
  ) |>
    filter(!grepl("output_manifest\\.csv$", file))
  write_csv(manifest, file.path(OUTPUT_DIRS$metadata, "output_manifest.csv"))

  invisible(results)
}
