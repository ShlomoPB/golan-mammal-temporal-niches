run_part6 <- function(data, part4) {
  message("Part 6: gazelle human-activity and reserve focus")
  focal <- data$focal
  captures <- data$captures

  gazelle_human <- part4$human_summary |>
    filter(species == "Gazella gazella")
  gazelle_curves <- part4$curves |>
    filter(species == "Gazella gazella")

  panel_a_overlap <- gazelle_curves |>
    group_by(hour) |>
    summarise(overlap_density = min(density), .groups = "drop")

  n_low <- unique(gazelle_curves$n_low)
  n_high <- unique(gazelle_curves$n_high)
  panel_a_labels <- c(
    "Low human activity (Q1)" = paste0(
      "Q1 low human activity (n=", format(n_low, big.mark = ","), ")"
    ),
    "High human activity (Q4)" = paste0(
      "Q4 high human activity (n=", format(n_high, big.mark = ","), ")"
    )
  )

  panel_a <- ggplot(gazelle_curves, aes(hour, density, color = group)) +
    daypart_layers() +
    geom_area(
      data = panel_a_overlap,
      aes(x = hour, y = overlap_density),
      inherit.aes = FALSE,
      fill = "grey45",
      alpha = 0.23
    ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = HUMAN_ACTIVITY_COLORS,
      labels = panel_a_labels,
      name = NULL
    ) +
    scale_x_continuous(
      breaks = c(0, 6, 12, 18, 24),
      labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
      limits = c(-0.25, 24.35),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.10, by = 0.02),
      labels = label_number(accuracy = 0.01),
      limits = c(0, 0.105),
      expand = c(0, 0)
    ) +
    labs(
      title = "(A) All land-use types",
      x = "Solar time",
      y = "Mountain gazelle temporal activity density"
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 19),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 19),
      axis.text = element_text(size = 17),
      legend.position = "inside",
      legend.position.inside = c(0.985, 0.985),
      legend.justification = c(1, 1),
      legend.background = element_rect(
        fill = alpha("white", 0.92),
        color = "grey55"
      ),
      legend.text = element_text(size = 16),
      plot.margin = margin(10, 18, 8, 8)
    )

  reserve_data <- captures |>
    mutate(
      reserve = case_when(
        camera_name %in% as.character(6:10) ~ "Temporally zoned reserve",
        camera_name %in% as.character(31:35) ~ "Open-access reserve",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(reserve), is.finite(time_rad))

  gazelle_reserve <- reserve_data |>
    filter(species == "Gazella gazella")
  reserve_counts <- gazelle_reserve |>
    count(reserve, name = "n_events")

  reserve_curves <- map_dfr(names(RESERVE_COLORS), function(group) {
    curve <- activity_curve(
      gazelle_reserve |> filter(reserve == group),
      weighted = FALSE
    )
    curve |> mutate(reserve = group)
  })
  human_curves <- map_dfr(names(RESERVE_COLORS), function(group) {
    curve <- activity_curve(
      reserve_data |> filter(reserve == group, species %in% HUMAN_ALIASES),
      weighted = FALSE
    )
    curve |> mutate(reserve = group)
  })

  reserve_overlap_fill <- reserve_curves |>
    group_by(hour) |>
    summarise(overlap_density = min(density), .groups = "drop")

  reserve_labels <- setNames(
    map_chr(names(RESERVE_COLORS), function(group) {
      n <- reserve_counts$n_events[reserve_counts$reserve == group]
      paste0(group, " (n=", format(n, big.mark = ","), ")")
    }),
    names(RESERVE_COLORS)
  )

  panel_b <- ggplot() +
    daypart_layers() +
    geom_area(
      data = reserve_overlap_fill,
      aes(hour, overlap_density),
      fill = "grey45",
      alpha = 0.23
    ) +
    geom_line(
      data = reserve_curves,
      aes(hour, density, color = reserve),
      linewidth = 1.2
    ) +
    geom_line(
      data = human_curves,
      aes(hour, density, color = reserve, linetype = "Human activity"),
      linewidth = 1.05,
      show.legend = TRUE
    ) +
    scale_color_manual(
      values = RESERVE_COLORS,
      labels = reserve_labels,
      limits = names(RESERVE_COLORS),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("Human activity" = "dotted"),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = c(0, 6, 12, 18, 24),
      labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
      limits = c(-0.25, 24.35),
      expand = c(0, 0)
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    labs(
      title = "(B) Nature reserves",
      x = "Solar time",
      y = "Mountain gazelle temporal activity density"
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 19),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 19),
      axis.text = element_text(size = 17),
      legend.position = c(0.015, 0.985),
      legend.justification = c(0, 1),
      legend.background = element_blank(),
      legend.box = "vertical",
      legend.box.background = element_rect(
        fill = alpha("white", 0.92),
        color = "grey55"
      ),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.margin = margin(5, 7, 5, 7),
      legend.text = element_text(size = 16),
      plot.margin = margin(10, 18, 8, 8)
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(linetype = "solid", linewidth = 1.2)
      ),
      linetype = guide_legend(
        order = 2,
        override.aes = list(color = "grey25", linewidth = 1.05)
      )
    )

  figure6 <- panel_a / panel_b + plot_layout(heights = c(1, 1))
  save_figure(
    figure6,
    file.path(
      OUTPUT_DIRS$main_figures,
      "Figure_6_gazelle_human_and_reserve_activity"
    ),
    width = 11.5,
    height = 10.2
  )

  set.seed(SEED + 600L)
  yehudiya <- gazelle_reserve |> filter(reserve == "Temporally zoned reserve")
  odem <- gazelle_reserve |> filter(reserve == "Open-access reserve")
  reserve_test <- adaptive_overlap_test(
    yehudiya$time_rad,
    odem$time_rad,
    yehudiya$sighting_quantity,
    odem$sighting_quantity,
    dhat4_reps = N_RANDOMIZATIONS_PART6,
    dhat1_reps = N_RANDOMIZATIONS_DHAT1,
    bootstrap_reps = N_BOOT_DHAT1
  ) |>
    mutate(
      comparison = "Temporally zoned reserve vs open-access reserve",
      n_temporally_zoned = nrow(yehudiya),
      n_open_access = nrow(odem),
      n_individuals_temporally_zoned = sum(yehudiya$sighting_quantity),
      n_individuals_open_access = sum(odem$sighting_quantity),
      .before = 1
    )
  write_csv(
    reserve_test,
    file.path(
      OUTPUT_DIRS$additional_results,
      "part6_gazelle_reserve_overlap_test.csv"
    )
  )

  gazelle_summary <- gazelle_human |>
    transmute(
      mean_noct_Q1,
      mean_noct_Q4,
      delta_noct_Q4_minus_Q1,
      delta_CI_low,
      delta_CI_high,
      n_gazelle_events_Q1 = n_low_Q1,
      n_gazelle_events_Q4 = n_high_Q4,
      n_gazelle_individuals_Q1 = n_individuals_low_Q1,
      n_gazelle_individuals_Q4 = n_individuals_high_Q4,
      n_cameras_Q1 = N_Q_EXTREME_CAMERAS,
      n_cameras_Q4 = N_Q_EXTREME_CAMERAS,
      overlap_estimator = estimator,
      overlap_delta = temporal_overlap_delta,
      overlap_p = kernel_test_p,
      nocturnality_p = nocturnality_test_p
    )
  write_csv(
    gazelle_summary,
    file.path(
      OUTPUT_DIRS$additional_results,
      "gazelle_human_activity_summary.csv"
    )
  )

  list(
    gazelle_summary = gazelle_summary,
    reserve_test = reserve_test,
    figure = figure6
  )
}
