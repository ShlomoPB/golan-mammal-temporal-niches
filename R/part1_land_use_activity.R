run_part1 <- function(data) {
  message("Part 1: land-use activity contrasts")
  focal <- data$focal
  land_pairs <- combn(LAND_USE_LEVELS, 2, simplify = FALSE)

  results <- map_dfr(SPECIES, function(sp) {
    message("  overlap tests: ", sp)
    sp_data <- focal |> filter(species == sp)
    map_dfr(land_pairs, function(pair) {
      group1 <- sp_data |> filter(land_use == pair[1])
      group2 <- sp_data |> filter(land_use == pair[2])
      test <- adaptive_overlap_test(
        group1$time_rad,
        group2$time_rad,
        group1$sighting_quantity,
        group2$sighting_quantity,
        dhat4_reps = N_RANDOMIZATIONS_PART1,
        dhat1_reps = N_RANDOMIZATIONS_DHAT1,
        bootstrap_reps = N_BOOT_DHAT1
      )
      test |>
        mutate(
          species = sp,
          group1 = pair[1],
          group2 = pair[2],
          comparison = paste(pair, collapse = " vs "),
          comparison_code = paste(LAND_USE_CODES[pair], collapse = "-"),
          .before = 1
        )
    })
  }) |>
    group_by(species) |>
    mutate(p_adj_holm = p.adjust(p_value, method = "holm")) |>
    ungroup() |>
    mutate(
      se = coalesce(bootstrap_se, null_se),
      se_source = if_else(
        estimator == "Dhat1",
        "within-group bootstrap SE",
        "compareCkern randomization-null SE"
      )
    )

  write_csv(
    results,
    file.path(OUTPUT_DIRS$additional_results, "part1_all_land_use_overlap_tests.csv")
  )

  table1 <- results |>
    transmute(
      Species = species,
      `Land-use pair` = comparison_code,
      `Group 1` = group1,
      `Group 2` = group2,
      `n events, group 1` = n1,
      `n events, group 2` = n2,
      Estimator = estimator,
      Overlap = delta,
      SE = se,
      `SE source` = se_source,
      p = p_value,
      p_adj = p_adj_holm
    )
  write_csv(table1, file.path(OUTPUT_DIRS$main_tables, "Table_1_land_use_activity_tests.csv"))

  curves <- map_dfr(SPECIES, function(sp) {
    map_dfr(LAND_USE_LEVELS, function(lu) {
      curve <- activity_curve(focal |> filter(species == sp, land_use == lu), weighted = TRUE)
      if (is.null(curve)) return(tibble())
      curve |> mutate(species = sp, land_use = lu)
    })
  })

  curve_areas <- curves |>
    arrange(species, land_use, hour) |>
    group_by(species, land_use) |>
    summarise(area = trapz(hour, density), .groups = "drop")
  write_csv(
    curve_areas,
    file.path(OUTPUT_DIRS$diagnostics, "part1_activity_curve_area_audit.csv")
  )

  selected <- results |>
    group_by(species) |>
    slice_min(delta, n = 1, with_ties = FALSE) |>
    ungroup()

  forced <- results |>
    filter(
      species == FIGURE3_FORCED_PAIR[["species"]],
      group1 == FIGURE3_FORCED_PAIR[["group1"]],
      group2 == FIGURE3_FORCED_PAIR[["group2"]]
    )
  selected <- selected |>
    filter(species != FIGURE3_FORCED_PAIR[["species"]]) |>
    bind_rows(forced) |>
    mutate(species = factor(species, levels = SPECIES)) |>
    arrange(species) |>
    mutate(
      selection_rule = if_else(
        as.character(species) == FIGURE3_FORCED_PAIR[["species"]],
        "Author-requested high-culling vs low-culling contrast",
        "Minimum estimated overlap"
      )
    )
  write_csv(
    selected,
    file.path(OUTPUT_DIRS$additional_results, "Figure_3_selected_land_use_pairs.csv")
  )

  plot_data <- map_dfr(seq_len(nrow(selected)), function(i) {
    row <- selected[i, ]
    curves |>
      filter(
        species == as.character(row$species),
        land_use %in% c(row$group1, row$group2)
      ) |>
      mutate(
        species_order = as.integer(row$species),
        panel = paste0(
          letters[i], ". <b><i>",
          as.character(row$species),
          "</i></b><br>",
          row$group1,
          " vs ",
          row$group2
        )
      )
  })

  overlap_fill <- plot_data |>
    group_by(species, species_order, panel, hour) |>
    summarise(
      overlap_density = min(density),
      n_curves = n_distinct(land_use),
      .groups = "drop"
    ) |>
    filter(n_curves == 2)

  plot_data <- plot_data |>
    mutate(
      panel = factor(panel, levels = unique(plot_data$panel)),
      land_use = factor(land_use, levels = LAND_USE_LEVELS)
    )
  overlap_fill$panel <- factor(overlap_fill$panel, levels = levels(plot_data$panel))

  figure3 <- ggplot(plot_data, aes(hour, density, color = land_use)) +
    daypart_layers() +
    geom_area(
      data = overlap_fill,
      aes(x = hour, y = overlap_density, group = panel),
      inherit.aes = FALSE,
      fill = "grey45",
      alpha = 0.23
    ) +
    geom_line(linewidth = 1.05) +
    facet_wrap(~panel, ncol = 3, scales = "free_y") +
    scale_color_manual(
      values = LAND_USE_COLORS,
      limits = LAND_USE_LEVELS,
      drop = FALSE,
      name = "Land-use category"
    ) +
    scale_x_continuous(
      breaks = c(0, 6, 12, 18, 24),
      labels = c("00:00", "06:00", "12:00", "18:00", "24:00"),
      limits = c(-0.25, 24.35),
      expand = c(0, 0)
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    labs(x = "Solar time", y = "Temporal activity density") +
    guides(
      color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 1.8)
      )
    ) +
    theme_bw(base_size = 16) +
    theme(
      strip.text = ggtext::element_markdown(size = 12.5, lineheight = 1.08),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.70, 0.19),
      legend.justification = c(0.5, 0.5),
      legend.direction = "horizontal",
      legend.title = element_text(face = "plain", size = 20),
      legend.text = element_text(size = 18),
      legend.key.width = grid::unit(2.2, "cm"),
      legend.spacing.x = grid::unit(0.55, "cm"),
      legend.background = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      plot.margin = margin(8, 16, 8, 8)
    )

  save_figure(
    figure3,
    file.path(OUTPUT_DIRS$main_figures, "Figure_3_land_use_activity"),
    width = 13.2,
    height = 9.3
  )

  list(tests = results, table1 = table1, curves = curves, selected = selected)
}
