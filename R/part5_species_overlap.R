run_part5 <- function(data) {
  message("Part 5: interspecific temporal overlap")
  focal <- data$focal
  species_pairs <- combn(SPECIES, 2, simplify = FALSE)

  set.seed(SEED + 500L)
  overlap_results <- map_dfr(LAND_USE_LEVELS, function(lu) {
    message("  bootstrap overlap CIs: ", lu)
    treatment_data <- focal |> filter(land_use == lu)
    map_dfr(species_pairs, function(pair) {
      group1 <- treatment_data |> filter(species == pair[1])
      group2 <- treatment_data |> filter(species == pair[2])
      estimator <- choose_delta_estimator(nrow(group1), nrow(group2))
      boot <- bootstrap_standard_delta(
        group1$time_rad,
        group2$time_rad,
        estimator,
        reps = N_BOOT_INTERSPECIES_OVERLAP
      )
      boot |>
        mutate(
          land_use = lu,
          species1 = pair[1],
          species2 = pair[2],
          n1 = nrow(group1),
          n2 = nrow(group2),
          n_individuals1 = sum(group1$sighting_quantity),
          n_individuals2 = sum(group2$sighting_quantity),
          estimator = estimator,
          .before = 1
        )
    })
  })

  write_csv(
    overlap_results,
    file.path(
      OUTPUT_DIRS$additional_results,
      "part5_pairwise_species_overlap_long.csv"
    )
  )

  land_use_table_order <- c(
    "High culling",
    "Low culling",
    "Minefields",
    "Nature reserves"
  )
  table_a4 <- overlap_results |>
    mutate(
      land_use = factor(land_use, levels = land_use_table_order),
      Species1 = unname(SPECIES_ABBREVIATED[species1]),
      Species2 = unname(SPECIES_ABBREVIATED[species2])
    ) |>
    select(
      Species1,
      Species2,
      land_use,
      Estimator = estimator,
      Overlap = delta,
      CI_low = ci_low,
      CI_high = ci_high,
      N1 = n1,
      N2 = n2
    ) |>
    pivot_wider(
      names_from = land_use,
      values_from = c(Estimator, Overlap, CI_low, CI_high, N1, N2),
      names_glue = "{land_use}_{.value}",
      names_vary = "slowest"
    )
  write_csv(
    table_a4,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_A4_pairwise_overlap_by_land_use.csv"
    )
  )

  model_data <- overlap_results |>
    mutate(
      species_a = pmin(species1, species2),
      species_b = pmax(species1, species2),
      pair_id = paste(species_a, species_b, sep = "~"),
      land_use = factor(
        land_use,
        levels = c("Minefields", "Nature reserves", "Low culling", "High culling")
      ),
      disturbance_rank = as.numeric(land_use) - 1,
      is_gazelle_pair =
        species1 == "Gazella gazella" | species2 == "Gazella gazella"
    )
  n_all <- nrow(model_data)
  model_data <- model_data |>
    mutate(overlap_beta = (delta * (n_all - 1) + 0.5) / n_all)

  extract_models <- function(df, subset_label) {
    rank_model <- glmmTMB(
      overlap_beta ~ disturbance_rank + (1 | pair_id),
      family = beta_family(link = "logit"),
      data = df
    )
    unranked_model <- glmmTMB(
      overlap_beta ~ land_use + (1 | pair_id),
      family = beta_family(link = "logit"),
      data = df
    )
    null_model <- glmmTMB(
      overlap_beta ~ 1 + (1 | pair_id),
      family = beta_family(link = "logit"),
      data = df
    )

    rank_coef <- summary(rank_model)$coefficients$cond["disturbance_rank", ]
    rank_ci <- suppressMessages(confint(rank_model, parm = "disturbance_rank"))
    rank_row <- tibble(
      subset = subset_label,
      section = "Rank slope (logit)",
      term = "disturbance_rank",
      estimate = unname(rank_coef["Estimate"]),
      SE = unname(rank_coef["Std. Error"]),
      statistic = unname(rank_coef["z value"]),
      p_value = unname(rank_coef["Pr(>|z|)"]),
      CI_low = unname(rank_ci[1, 1]),
      CI_high = unname(rank_ci[1, 2])
    )

    emm_rank <- as.data.frame(summary(
      emmeans(
        rank_model,
        ~disturbance_rank,
        at = list(disturbance_rank = 0:3),
        type = "response"
      )
    ))
    rank_emm_rows <- tibble(
      subset = subset_label,
      section = "Rank-model EMM",
      term = c("Minefields", "Nature reserves", "Low culling", "High culling"),
      estimate = emm_rank$response,
      SE = emm_rank$SE,
      statistic = NA_real_,
      p_value = NA_real_,
      CI_low = emm_rank$asymp.LCL,
      CI_high = emm_rank$asymp.UCL
    )

    emm_unranked <- emmeans(unranked_model, ~land_use, type = "response")
    emm_unranked_df <- as.data.frame(summary(emm_unranked))
    unranked_emm_rows <- tibble(
      subset = subset_label,
      section = "Unranked EMM",
      term = as.character(emm_unranked_df$land_use),
      estimate = emm_unranked_df$response,
      SE = emm_unranked_df$SE,
      statistic = NA_real_,
      p_value = NA_real_,
      CI_low = emm_unranked_df$asymp.LCL,
      CI_high = emm_unranked_df$asymp.UCL
    )

    pair_df <- as.data.frame(summary(
      pairs(emm_unranked, adjust = "holm"),
      infer = c(TRUE, TRUE)
    ))
    estimate_column <- intersect(
      c("odds.ratio", "ratio", "estimate"),
      names(pair_df)
    )[1]
    low_column <- intersect(c("asymp.LCL", "lower.CL"), names(pair_df))[1]
    high_column <- intersect(c("asymp.UCL", "upper.CL"), names(pair_df))[1]
    contrast_rows <- tibble(
      subset = subset_label,
      section = "Unranked contrasts",
      term = pair_df$contrast,
      estimate = pair_df[[estimate_column]],
      SE = pair_df$SE,
      statistic = pair_df$z.ratio,
      p_value = pair_df$p.value,
      CI_low = pair_df[[low_column]],
      CI_high = pair_df[[high_column]]
    )

    wald <- car::Anova(unranked_model, type = "III")
    wald_row <- tibble(
      subset = subset_label,
      section = "Global test",
      term = "Wald chi-square: land use",
      estimate = unname(wald["land_use", "Chisq"]),
      SE = NA_real_,
      statistic = unname(wald["land_use", "Df"]),
      p_value = unname(wald["land_use", "Pr(>Chisq)"]),
      CI_low = NA_real_,
      CI_high = NA_real_
    )

    lrt <- anova(null_model, unranked_model)
    lrt_row <- tibble(
      subset = subset_label,
      section = "Global test",
      term = "Likelihood-ratio test vs null",
      estimate = unname(lrt$Chisq[2]),
      SE = NA_real_,
      statistic = unname(lrt$`Chi Df`[2]),
      p_value = unname(lrt$`Pr(>Chisq)`[2]),
      CI_low = NA_real_,
      CI_high = NA_real_
    )

    list(
      table = bind_rows(
        rank_row,
        rank_emm_rows,
        unranked_emm_rows,
        contrast_rows,
        wald_row,
        lrt_row
      ),
      rank_model = rank_model,
      unranked_model = unranked_model,
      null_model = null_model
    )
  }

  model_non_gazelle <- extract_models(
    model_data |> filter(!is_gazelle_pair),
    "Non-gazelle pairs"
  )
  model_gazelle <- extract_models(
    model_data |> filter(is_gazelle_pair),
    "Gazelle vs others"
  )
  model_diagnostics <- list(
    "Non-gazelle pairs: rank" = model_non_gazelle$rank_model,
    "Non-gazelle pairs: unranked" = model_non_gazelle$unranked_model,
    "Non-gazelle pairs: null" = model_non_gazelle$null_model,
    "Gazelle vs others: rank" = model_gazelle$rank_model,
    "Gazelle vs others: unranked" = model_gazelle$unranked_model,
    "Gazelle vs others: null" = model_gazelle$null_model
  ) |>
    imap_dfr(function(fitted_model, model_name) {
      tibble(
        model = model_name,
        optimizer_convergence_code = fitted_model$fit$convergence,
        optimizer_message = fitted_model$fit$message,
        positive_definite_hessian = isTRUE(fitted_model$sdr$pdHess)
      )
    })
  write_csv(
    model_diagnostics,
    file.path(
      OUTPUT_DIRS$diagnostics,
      "part5_beta_model_convergence.csv"
    )
  )

  table_a5 <- bind_rows(model_non_gazelle$table, model_gazelle$table)
  write_csv(
    table_a5,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_A5_pairwise_overlap_beta_models.csv"
    )
  )

  bin_labels <- c("<0.4", "0.4-0.6", "0.6-0.8", "0.8-0.9", "0.9-1.0")
  bin_colors <- c(
    "<0.4" = "#FFFFB2",
    "0.4-0.6" = "#FECC5C",
    "0.6-0.8" = "#FD8D3C",
    "0.8-0.9" = "#F03B20",
    "0.9-1.0" = "#BD0026"
  )
  plot_data <- overlap_results |>
    mutate(
      x = match(species1, SPECIES),
      y = match(species2, SPECIES),
      overlap_bin = cut(
        delta,
        breaks = c(-Inf, 0.4, 0.6, 0.8, 0.9, Inf),
        labels = bin_labels,
        right = FALSE
      )
    )

  panel_letters <- setNames(c("A", "B", "C", "D"), LAND_USE_LEVELS)
  make_triangle <- function(lu, include_legend = TRUE) {
    panel <- ggplot(plot_data |> filter(land_use == lu), aes(x, y, fill = overlap_bin)) +
      geom_tile(color = "white", linewidth = 1.0, width = 0.98, height = 0.98) +
      scale_fill_manual(
        values = bin_colors,
        limits = bin_labels,
        drop = FALSE,
        name = expression(Delta)
      ) +
      scale_x_continuous(
        breaks = seq_along(SPECIES),
        labels = unname(SPECIES_ABBREVIATED[SPECIES]),
        limits = c(0.5, length(SPECIES) + 0.5),
        expand = c(0, 0)
      ) +
      scale_y_reverse(
        breaks = seq_along(SPECIES),
        labels = unname(SPECIES_ABBREVIATED[SPECIES]),
        limits = c(length(SPECIES) + 0.5, 0.5),
        expand = c(0, 0)
      ) +
      coord_fixed(clip = "off") +
      labs(title = paste0("(", panel_letters[lu], ") ", lu), x = NULL, y = NULL) +
      theme_void(base_size = 15) +
      theme(
        plot.title = element_text(face = "bold", size = 19, hjust = 0.5),
        axis.text.x = element_text(
          face = "italic",
          color = "black",
          size = 13,
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        axis.text.y = element_text(face = "italic", color = "black", size = 13),
        legend.position = if (include_legend) "bottom" else "none",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.key.height = grid::unit(0.55, "cm"),
        legend.key.width = grid::unit(0.65, "cm"),
        plot.margin = margin(12, 18, 12, 18)
      ) +
      guides(fill = guide_legend(nrow = 1, title.position = "left"))
    panel
  }

  triangle_plots <- setNames(
    map(LAND_USE_LEVELS, ~make_triangle(.x, include_legend = TRUE)),
    LAND_USE_LEVELS
  )
  walk2(triangle_plots, LAND_USE_LEVELS, function(plot, lu) {
    file_stub <- gsub(" ", "_", lu)
    save_figure(
      plot,
      file.path(
        OUTPUT_DIRS$main_figures,
        paste0("Figure_5_", file_stub)
      ),
      width = 8.1,
      height = 8.3
    )
  })

  combined_panels <- wrap_plots(
    map(LAND_USE_LEVELS, ~make_triangle(.x, include_legend = FALSE)),
    ncol = 2
  )
  legend_plot <- ggplot(
    tibble(
      x = seq_along(bin_labels),
      y = 1,
      overlap_bin = factor(bin_labels, levels = bin_labels)
    ),
    aes(x, y, fill = overlap_bin)
  ) +
    geom_point(shape = 22, size = 0, alpha = 0) +
    scale_fill_manual(
      values = bin_colors,
      limits = bin_labels,
      drop = FALSE,
      name = expression(Delta)
    ) +
    guides(fill = guide_legend(
      nrow = 1,
      title.position = "left",
      override.aes = list(size = 6, alpha = 1, shape = 22)
    )) +
    theme_void(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.key.height = grid::unit(0.55, "cm"),
      legend.key.width = grid::unit(0.65, "cm")
    )
  combined <- combined_panels / legend_plot +
    plot_layout(heights = c(1, 0.07))
  save_figure(
    combined,
    file.path(OUTPUT_DIRS$main_figures, "Figure_5_pairwise_overlap_matrices"),
    width = 13.2,
    height = 12.4
  )

  list(
    overlaps = overlap_results,
    table_a4 = table_a4,
    table_a5 = table_a5,
    model_diagnostics = model_diagnostics,
    models = list(
      non_gazelle = model_non_gazelle,
      gazelle = model_gazelle
    )
  )
}
