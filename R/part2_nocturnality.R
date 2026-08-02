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
  forest_data <- slopes |>
    mutate(
      species = factor(species, levels = rev(c(SPECIES, "Community"))),
      species_label = if_else(
        as.character(species) == "Community",
        "<b>Community</b>",
        paste0("<i>", as.character(species), "</i>")
      ),
      species_label = factor(
        species_label,
        levels = c("<b>Community</b>", paste0("<i>", rev(SPECIES), "</i>"))
      ),
      group = case_when(
        as.character(species) == "Community" ~ "Community",
        as.character(species) %in% culled_species ~ "Culled species",
        TRUE ~ "Other species"
      )
    )
  forest_colors <- c(
    "Community" = "black",
    "Culled species" = "#B85A68",
    "Other species" = "grey35"
  )

  figure4 <- ggplot(forest_data, aes(estimate, species_label, color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = 0,
      linewidth = 0.8
    ) +
    geom_point(size = 3.2) +
    scale_color_manual(values = forest_colors, guide = "none") +
    labs(
      x = "Change in relative nocturnality per land-use step (h)",
      y = NULL
    ) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.y = ggtext::element_markdown(),
      plot.margin = margin(8, 14, 8, 8)
    )
  save_figure(
    figure4,
    file.path(OUTPUT_DIRS$main_figures, "Figure_4_nocturnality_gradient"),
    width = 8.2,
    height = 6.4
  )

  figure_a1_data <- species_land_use |>
    mutate(
      species_label = paste0("<b><i>", as.character(species), "</i></b>"),
      species_label = factor(
        species_label,
        levels = paste0("<b><i>", SPECIES, "</i></b>")
      )
    )

  figure_a1 <- ggplot(
    figure_a1_data,
    aes(land_use, mean, group = species)
  ) +
    geom_hline(yintercept = seq(0, 12, 2), color = "grey92", linewidth = 0.35) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.08, linewidth = 0.55) +
    geom_line(color = "grey35", linewidth = 0.7) +
    geom_point(size = 2.2, color = "#3F7FB5") +
    facet_wrap(~species_label, ncol = 3) +
    scale_x_discrete(labels = c(
      "High culling" = "High\nculling",
      "Low culling" = "Low\nculling",
      "Nature reserves" = "Nature\nreserves",
      "Minefields" = "Minefields"
    )) +
    coord_cartesian(ylim = c(4.5, 10)) +
    labs(x = "Land-use category", y = "Relative nocturnality (h from solar noon)") +
    theme_bw(base_size = 14) +
    theme(
      strip.text = ggtext::element_markdown(size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 14)
    )
  save_figure(
    figure_a1,
    file.path(OUTPUT_DIRS$supplementary_figures, "Figure_A1_nocturnality_by_land_use"),
    width = 11.5,
    height = 8.2
  )

  list(
    species_land_use = species_land_use,
    community_land_use = community_land_use,
    slopes = slopes
  )
}
