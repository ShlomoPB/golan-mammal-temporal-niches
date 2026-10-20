run_part4 <- function(data) {
  message("Part 4: species responses to camera-level pressures")
  focal <- data$focal
  site_covs <- data$site_covs

  set.seed(SEED + 400L)
  tests <- pmap_dfr(PART4_SPECIES_PRESSURES, function(species, pressure) {
    covariate <- unname(PART4_PRESSURES[pressure])
    groups <- select_extreme_cameras(site_covs, covariate)
    species_data <- focal |> filter(.data$species == .env$species)
    low <- species_data |> filter(camera_name %in% groups$low)
    high <- species_data |> filter(camera_name %in% groups$high)

    adaptive_overlap_test(
      low$time_rad,
      high$time_rad,
      low$sighting_quantity,
      high$sighting_quantity,
      dhat4_reps = N_RANDOMIZATIONS_PART4,
      dhat1_reps = N_RANDOMIZATIONS_DHAT1,
      bootstrap_reps = N_BOOT_DHAT1
    ) |>
      mutate(
        species = species,
        pressure = pressure,
        covariate = covariate,
        n_cameras_low = length(groups$low),
        n_cameras_high = length(groups$high),
        .before = 1
      )
  }) |>
    group_by(species) |>
    mutate(p_adj_holm = p.adjust(p_value, method = "holm")) |>
    ungroup()

  write_csv(
    tests,
    file.path(OUTPUT_DIRS$additional_results, "part4_all_pressure_overlap_tests.csv")
  )

  human_groups <- select_extreme_cameras(site_covs, "human.abundance")
  write_csv(
    human_groups$assignments |>
      left_join(
        site_covs |> select(camera_name, land_use),
        by = "camera_name"
      ),
    file.path(
      OUTPUT_DIRS$diagnostics,
      "part4_human_activity_Q1_Q4_camera_assignments.csv"
    )
  )

  human_summary <- map_dfr(SPECIES, function(sp) {
    species_data <- focal |> filter(species == sp)
    low <- species_data |> filter(camera_name %in% human_groups$low)
    high <- species_data |> filter(camera_name %in% human_groups$high)

    low_mean <- weighted_mean_bootstrap(low$nocturnality_h, low$sighting_quantity)
    high_mean <- weighted_mean_bootstrap(high$nocturnality_h, high$sighting_quantity)
    difference <- weighted_difference_bootstrap(
      low$nocturnality_h,
      low$sighting_quantity,
      high$nocturnality_h,
      high$sighting_quantity
    )

    lm_data <- bind_rows(
      low |> mutate(human_group = 0),
      high |> mutate(human_group = 1)
    )
    noct_fit <- lm(
      nocturnality_h ~ human_group,
      data = lm_data,
      weights = sighting_quantity
    )
    noct_coef <- summary(noct_fit)$coefficients["human_group", ]

    overlap_row <- tests |>
      filter(species == sp, pressure == "human activity")

    tibble(
      species = sp,
      n_high_Q4 = nrow(high),
      n_low_Q1 = nrow(low),
      n_individuals_high_Q4 = sum(high$sighting_quantity),
      n_individuals_low_Q1 = sum(low$sighting_quantity),
      estimator = overlap_row$estimator,
      temporal_overlap_delta = overlap_row$delta,
      kernel_test_p = overlap_row$p_value,
      mean_noct_Q1 = low_mean$mean,
      mean_noct_Q4 = high_mean$mean,
      delta_noct_Q4_minus_Q1 = difference$delta,
      delta_CI_low = difference$delta_ci_low,
      delta_CI_high = difference$delta_ci_high,
      nocturnality_test_p = unname(noct_coef["Pr(>|t|)"])
    )
  })

  table_a1 <- human_summary |>
    transmute(
      Species = species,
      `n, High human activity (Q4)` = n_high_Q4,
      `n, Low human activity (Q1)` = n_low_Q1,
      Estimator = estimator,
      `Temporal overlap (Delta)` = temporal_overlap_delta,
      `Kernel test p` = kernel_test_p,
      `Mean nocturnality Q1 (h)` = mean_noct_Q1,
      `Mean nocturnality Q4 (h)` = mean_noct_Q4,
      `Delta nocturnality, Q4-Q1 (h)` = delta_noct_Q4_minus_Q1,
      `Delta 95% CI low` = delta_CI_low,
      `Delta 95% CI high` = delta_CI_high,
      `Nocturnality test p` = nocturnality_test_p
    )
  write_csv(
    table_a1,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_A1_species_responses_to_human_activity.csv"
    )
  )

  curves <- map_dfr(SPECIES, function(sp) {
    species_data <- focal |> filter(species == sp)
    low <- species_data |> filter(camera_name %in% human_groups$low)
    high <- species_data |> filter(camera_name %in% human_groups$high)
    low_curve <- activity_curve(low, weighted = TRUE)
    high_curve <- activity_curve(high, weighted = TRUE)
    bind_rows(
      low_curve |> mutate(group = "Low human activity (Q1)"),
      high_curve |> mutate(group = "High human activity (Q4)")
    ) |>
      mutate(
        species = sp,
        n_low = nrow(low),
        n_high = nrow(high)
      )
  }) |>
    mutate(
      species = factor(species, levels = SPECIES),
      group = factor(group, levels = names(HUMAN_ACTIVITY_COLORS)),
      panel = paste0(
        unname(SPECIES_COMMON[as.character(species)]),
        " (n = ",
        n_low,
        " / ",
        n_high,
        ")"
      )
    )
  panel_levels <- curves |>
    distinct(species, panel) |>
    arrange(species) |>
    pull(panel)
  curves <- curves |>
    mutate(panel = factor(panel, levels = panel_levels))

  overlap_fill <- curves |>
    group_by(species, panel, hour) |>
    summarise(overlap_density = min(density), .groups = "drop")

  significance_stars <- function(p) {
    case_when(
      is.finite(p) & p < 0.01 ~ "**",
      is.finite(p) & p < 0.05 ~ "*",
      TRUE ~ ""
    )
  }
  annotation_data <- human_summary |>
    transmute(
      species = factor(species, levels = SPECIES),
      annotation = paste0(
        "\u0394",
        if_else(estimator == "Dhat1", "\u2081", "\u2084"),
        " = ",
        sprintf("%.2f", temporal_overlap_delta),
        significance_stars(kernel_test_p),
        "\nNoct. = ",
        sprintf("%+.2f", delta_noct_Q4_minus_Q1),
        " h",
        significance_stars(nocturnality_test_p)
      )
    ) |>
    left_join(
      curves |>
        group_by(species, panel) |>
        summarise(curve_max = max(density), .groups = "drop"),
      by = "species"
    ) |>
    mutate(
      panel = factor(panel, levels = panel_levels),
      hour = 12,
      label_y = curve_max * 1.18
    )

  figure_a2 <- ggplot(curves, aes(hour, density, color = group)) +
    daypart_layers() +
    geom_area(
      data = overlap_fill,
      aes(x = hour, y = overlap_density, group = panel),
      inherit.aes = FALSE,
      fill = "grey45",
      alpha = 0.23
    ) +
    geom_line(linewidth = 1.0) +
    geom_label(
      data = annotation_data,
      aes(x = hour, y = label_y, label = annotation),
      inherit.aes = FALSE,
      color = "black",
      fill = alpha("white", 0.92),
      linewidth = 0.3,
      label.padding = grid::unit(0.19, "lines"),
      lineheight = 1,
      size = 4
    ) +
    facet_wrap(~panel, ncol = 3, scales = "free_y") +
    scale_color_manual(values = HUMAN_ACTIVITY_COLORS, name = NULL) +
    scale_x_continuous(
      breaks = c(0, 6, 12, 18, 24),
      labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
      limits = c(-0.25, 24.35),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.01, 0.08))
    ) +
    labs(
      x = "Solar time",
      y = "Temporal activity density"
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 1.6)
      )
    ) +
    theme_bw(base_size = 15) +
    theme(
      strip.text = element_text(size = 12.5, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.70, 0.18),
      legend.justification = c(0.5, 0.5),
      legend.direction = "vertical",
      legend.text = element_text(size = 16),
      legend.key.width = grid::unit(1.8, "cm"),
      legend.spacing.x = grid::unit(0.45, "cm"),
      legend.background = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 11),
      plot.margin = margin(8, 16, 8, 8)
    )
  save_figure(
    figure_a2,
    file.path(
      OUTPUT_DIRS$supplementary_figures,
      "Figure_A2_activity_by_human_activity"
    ),
    width = 13,
    height = 9.2
  )

  curve_audit <- curves |>
    group_by(species, group) |>
    summarise(area = trapz(hour, density), .groups = "drop")
  write_csv(
    curve_audit,
    file.path(OUTPUT_DIRS$diagnostics, "part4_activity_curve_area_audit.csv")
  )

  list(
    tests = tests,
    human_summary = human_summary,
    table_a1 = table_a1,
    curves = curves,
    human_groups = human_groups
  )
}
