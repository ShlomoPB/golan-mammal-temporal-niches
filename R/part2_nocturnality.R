run_part2 <- function(data) {
  message("Part 2: nocturnality across the land-use gradient")
  focal <- data$focal |>
    mutate(
      land_use_index = match(as.character(land_use), LAND_USE_GRADIENT)
    )

  set.seed(SEED + 200L)
  species_land_use <- focal |>
    group_by(species, land_use) |>
    group_modify(~{
      summary <- weighted_mean_bootstrap(
        .x$nocturnality_h,
        .x$sighting_quantity
      )
      summary |>
        mutate(
          n_events = nrow(.x),
          n_individuals = sum(.x$sighting_quantity)
        )
    }) |>
    ungroup() |>
    mutate(
      species = factor(species, levels = SPECIES),
      land_use = factor(land_use, levels = LAND_USE_LEVELS)
    ) |>
    arrange(species, land_use)

  table_a3 <- species_land_use |>
    transmute(
      Species = as.character(species),
      land_use = as.character(land_use),
      noct_h_mean = mean,
      CI_lo = ci_low,
      CI_hi = ci_high,
      n_events,
      n_indiv = n_individuals
    )
  write_csv(
    table_a3,
    file.path(OUTPUT_DIRS$supplementary_tables, "Table_A3_nocturnality_by_species_land_use.csv")
  )

  community_land_use <- focal |>
    group_by(land_use) |>
    group_modify(~{
      weighted_mean_bootstrap(.x$nocturnality_h, .x$sighting_quantity) |>
        mutate(n_events = nrow(.x), n_individuals = sum(.x$sighting_quantity))
    }) |>
    ungroup()
  write_csv(
    community_land_use,
    file.path(OUTPUT_DIRS$additional_results, "part2_community_nocturnality_by_land_use.csv")
  )

  fit_slope <- function(df, label) {
    fit <- lm(
      nocturnality_h ~ land_use_index,
      data = df,
      weights = sighting_quantity
    )
    coef_row <- summary(fit)$coefficients["land_use_index", ]
    ci <- confint(fit, "land_use_index", level = 0.95)
    tibble(
      species = label,
      estimate = unname(coef_row["Estimate"]),
      se = unname(coef_row["Std. Error"]),
      statistic = unname(coef_row["t value"]),
      p_value = unname(coef_row["Pr(>|t|)"]),
      ci_low = unname(ci[1]),
      ci_high = unname(ci[2]),
      n_events = nrow(df),
      n_individuals = sum(df$sighting_quantity)
    )
  }

  slopes <- bind_rows(
    map_dfr(SPECIES, function(sp) {
      fit_slope(focal |> filter(species == sp), sp)
    }),
    fit_slope(focal, "Community")
  )
  write_csv(
    slopes,
    file.path(OUTPUT_DIRS$additional_results, "part2_nocturnality_gradient_slopes.csv")
  )

  culled_species <- c("Canis lupus", "Canis aureus", "Sus scrofa")
  figure4_species_order <- c(
    "Community",
    "Sus scrofa",
    "Canis aureus",
    "Canis lupus",
    "Meles meles",
    "Hystrix indica",
    "Vulpes vulpes",
    "Gazella gazella"
  )
  figure4_labels <- c(
    "Community" = "Community (pooled)",
    SPECIES_COMMON
  )
  forest_data <- slopes |>
    mutate(
      species = factor(species, levels = figure4_species_order),
      species_label = unname(figure4_labels[as.character(species)]),
      species_label = factor(
        species_label,
        levels = rev(unname(figure4_labels[figure4_species_order]))
      ),
      group = case_when(
        as.character(species) == "Community" ~ "Community",
        as.character(species) %in% culled_species ~ "Culled species",
        TRUE ~ "Other species"
      ),
      group = factor(
        group,
        levels = c("Community", "Culled species", "Other species")
      )
    )
  forest_colors <- c(
    "Community" = "black",
    "Culled species" = "#3B82F6",
    "Other species" = "#9CA3AF"
  )

  figure4 <- ggplot(forest_data, aes(estimate, species_label, color = group)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.7
    ) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = 0,
      linewidth = 0.8
    ) +
    geom_point(size = 3.6) +
    scale_color_manual(values = forest_colors, name = NULL) +
    scale_x_continuous(
      breaks = seq(-0.75, 0.25, by = 0.25),
      labels = label_number(accuracy = 0.01),
      limits = c(-0.765, 0.27),
      expand = expansion(mult = 0)
    ) +
    labs(
      x = "Slope (hours change in nocturnality per land-use step)",
      y = NULL
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 0.8, size = 3.6)
      )
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text = element_text(color = "grey30"),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 17),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.55),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey45", linewidth = 0.65),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 14),
      legend.margin = margin(t = 5, r = 0, b = 0, l = 0),
      plot.margin = margin(8, 14, 4, 8)
    )
  save_figure(
    figure4,
    file.path(OUTPUT_DIRS$main_figures, "Figure_4_nocturnality_gradient"),
    width = 8.2,
    height = 6.4
  )

  # Figure A.1 presents effect sizes directly. Independent detection events
  # are sampled uniformly within each group, and group size is used once as
  # the analysis weight in each bootstrap weighted mean. Contrast intervals
  # come from replicate contrasts, not subtraction of group-mean CI limits.
  set.seed(SEED + 201L)
  figure_a1_data <- map_dfr(SPECIES, function(sp) {
    species_data <- focal |>
      filter(as.character(species) == sp)
    minefield_data <- species_data |>
      filter(as.character(land_use) == "Minefields")

    map_dfr(LAND_USE_LEVELS, function(lu) {
      land_use_data <- species_data |>
        filter(as.character(land_use) == lu)
      mean_nocturnality <- weighted_mean(
        land_use_data$nocturnality_h,
        land_use_data$sighting_quantity
      )

      if (lu == "Minefields") {
        contrast <- tibble(
          delta = 0,
          delta_ci_low = NA_real_,
          delta_ci_high = NA_real_
        )
      } else {
        contrast <- weighted_difference_bootstrap(
          minefield_data$nocturnality_h,
          minefield_data$sighting_quantity,
          land_use_data$nocturnality_h,
          land_use_data$sighting_quantity
        )
      }

      tibble(
        species_source = sp,
        species = unname(SPECIES_COMMON[[sp]]),
        land_use = lu,
        mean_nocturnality = mean_nocturnality,
        contrast_vs_minefields = contrast$delta,
        ci_low = contrast$delta_ci_low,
        ci_high = contrast$delta_ci_high,
        n_detections = nrow(land_use_data),
        sum_group_size_weights = sum(land_use_data$sighting_quantity)
      )
    })
  }) |>
    mutate(
      land_use = factor(land_use, levels = LAND_USE_LEVELS),
      species_label = paste0("<b>", species, "</b>"),
      species_label = factor(
        species_label,
        levels = paste0(
          "<b>",
          unname(SPECIES_COMMON[SPECIES]),
          "</b>"
        )
      )
    )

  mean_validation <- figure_a1_data |>
    left_join(
      species_land_use |>
        transmute(
          species_source = as.character(species),
          land_use = as.character(land_use),
          corrected_mean = mean,
          corrected_n_detections = n_events,
          corrected_sum_group_size_weights = n_individuals
        ),
      by = c("species_source", "land_use")
    )
  if (any(abs(
    mean_validation$mean_nocturnality - mean_validation$corrected_mean
  ) > 1e-12)) {
    stop(
      "Figure A.1 point estimates do not match the final corrected land-use means; ",
      "the figure was not saved."
    )
  }
  if (any(
    mean_validation$n_detections != mean_validation$corrected_n_detections |
      mean_validation$sum_group_size_weights !=
        mean_validation$corrected_sum_group_size_weights
  )) {
    stop(
      "Figure A.1 sample sizes do not match Table A.3; the figure was not saved."
    )
  }

  expected_rounded_contrasts <- tribble(
    ~species_source, ~land_use, ~expected_contrast,
    "Canis aureus", "High culling", 0.97,
    "Canis aureus", "Low culling", 0.51,
    "Canis aureus", "Nature reserves", 0.28,
    "Canis aureus", "Minefields", 0,
    "Canis lupus", "High culling", 0.54,
    "Canis lupus", "Low culling", 0.41,
    "Canis lupus", "Nature reserves", -0.04,
    "Canis lupus", "Minefields", 0,
    "Gazella gazella", "High culling", -0.31,
    "Gazella gazella", "Low culling", -0.26,
    "Gazella gazella", "Nature reserves", -0.05,
    "Gazella gazella", "Minefields", 0,
    "Hystrix indica", "High culling", 0.24,
    "Hystrix indica", "Low culling", 0.44,
    "Hystrix indica", "Nature reserves", 0.25,
    "Hystrix indica", "Minefields", 0,
    "Meles meles", "High culling", 0.36,
    "Meles meles", "Low culling", 0.11,
    "Meles meles", "Nature reserves", 0.56,
    "Meles meles", "Minefields", 0,
    "Sus scrofa", "High culling", 1.75,
    "Sus scrofa", "Low culling", 1.42,
    "Sus scrofa", "Nature reserves", 0.57,
    "Sus scrofa", "Minefields", 0,
    "Vulpes vulpes", "High culling", 0.15,
    "Vulpes vulpes", "Low culling", -0.03,
    "Vulpes vulpes", "Nature reserves", 0.25,
    "Vulpes vulpes", "Minefields", 0
  )
  contrast_discrepancies <- figure_a1_data |>
    mutate(land_use = as.character(land_use)) |>
    inner_join(
      expected_rounded_contrasts,
      by = c("species_source", "land_use")
    ) |>
    filter(abs(contrast_vs_minefields - expected_contrast) > 0.015)
  if (nrow(contrast_discrepancies)) {
    print(contrast_discrepancies, n = Inf)
    stop(
      "Figure A.1 contrasts differ from the final manuscript checks; ",
      "the figure was not saved."
    )
  }

  figure_a1_validation <- figure_a1_data |>
    transmute(
      species,
      land_use_category = as.character(land_use),
      weighted_mean_nocturnality = mean_nocturnality,
      contrast_vs_minefields,
      lower_95_CI = ci_low,
      upper_95_CI = ci_high,
      number_of_detections = n_detections,
      sum_of_group_size_weights = sum_group_size_weights
    )
  write_csv(
    figure_a1_validation,
    file.path(
      OUTPUT_DIRS$diagnostics,
      "Figure_A1_nocturnality_contrasts_validation.csv"
    )
  )
  estimates_outside_ci <- figure_a1_validation |>
    filter(
      land_use_category != "Minefields",
      contrast_vs_minefields < lower_95_CI |
        contrast_vs_minefields > upper_95_CI
    )
  if (nrow(estimates_outside_ci)) {
    print(estimates_outside_ci, n = Inf)
    stop(
      "Figure A.1 bootstrap percentile CI excludes the observed contrast for: ",
      paste(
        paste0(
          estimates_outside_ci$species,
          " (",
          estimates_outside_ci$land_use_category,
          ")"
        ),
        collapse = "; "
      ),
      ". The figure was not saved; review the estimator and data."
    )
  }
  message(
    "  Figure A.1 validation passed: all means match the corrected Part 2 ",
    "values, all contrasts match the manuscript checks, and every ",
    "non-reference point lies within its bootstrap CI."
  )

  figure_a1 <- ggplot(
    figure_a1_data,
    aes(land_use, contrast_vs_minefields, group = species)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.55
    ) +
    geom_errorbar(
      data = figure_a1_data |>
        filter(as.character(land_use) != "Minefields"),
      aes(ymin = ci_low, ymax = ci_high),
      width = 0.08,
      linewidth = 0.65
    ) +
    geom_line(color = "grey35", linewidth = 0.7) +
    geom_point(size = 2.5, color = "#3F7FB5") +
    facet_wrap(~species_label, ncol = 3) +
    scale_x_discrete(labels = c(
      "High culling" = "High\nculling",
      "Low culling" = "Low\nculling",
      "Nature reserves" = "Nature\nreserves",
      "Minefields" = "Minefields"
    )) +
    labs(
      x = "Land-use category",
      y = "Change in nocturnality relative to minefields (h)",
      subtitle = paste(
        "Points show mean change in nocturnality relative to minefields;",
        "error bars indicate 95% bootstrap confidence intervals."
      )
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.text = ggtext::element_markdown(size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.subtitle = element_text(size = 12.5, margin = margin(b = 8))
    )

  figure_a1_caption <- paste(
    "Figure A.1. Within-species change in nocturnality across land-use types.",
    "Points show mean change in nocturnality relative to minefields, measured",
    "as absolute hours from solar noon. Positive values indicate greater",
    "nocturnality than in minefields and negative values indicate less",
    "nocturnality. Error bars indicate 95% bootstrap confidence intervals.",
    "Minefields are the reference category and are therefore fixed at zero."
  )
  writeLines(
    figure_a1_caption,
    file.path(
      OUTPUT_DIRS$supplementary_figures,
      "Figure_A1_caption.txt"
    ),
    useBytes = TRUE
  )
  save_figure(
    figure_a1,
    file.path(OUTPUT_DIRS$supplementary_figures, "Figure_A1_nocturnality_by_land_use"),
    width = 11.5,
    height = 8.2
  )

  list(
    species_land_use = species_land_use,
    figure_a1_contrasts = figure_a1_data,
    community_land_use = community_land_use,
    slopes = slopes
  )
}
