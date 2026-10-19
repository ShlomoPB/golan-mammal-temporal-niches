suppressPackageStartupMessages({
  library(activity)
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(glmmTMB)
  library(janitor)
  library(overlap)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
  library(tibble)
  library(tidyr)
})

assert_file <- function(path) {
  if (!file.exists(path)) stop("Required input file is missing: ", path)
  invisible(path)
}

trapz <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2, na.rm = TRUE)
}

normalize_curve <- function(x, y) {
  area <- trapz(x, y)
  if (!is.finite(area) || area <= 0) stop("Activity curve has non-positive area.")
  y / area
}

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_mean_bootstrap <- function(x, w, reps = N_BOOT_NOCTURNALITY) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- as.numeric(x[keep])
  w <- as.numeric(w[keep])
  if (!length(x)) {
    return(tibble(mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  estimate <- weighted_mean(x, w)
  n <- length(x)
  draws <- replicate(reps, {
    # Independent detection events are the resampling unit. Group size remains
    # a weight in the statistic and is therefore not also used as a sampling
    # probability.
    id <- sample.int(n, n, replace = TRUE)
    weighted_mean(x[id], w[id])
  })

  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble(mean = estimate, ci_low = ci[1], ci_high = ci[2])
}

weighted_difference_bootstrap <- function(
    x_low, w_low, x_high, w_high, reps = N_BOOT_NOCTURNALITY) {
  keep_low <- is.finite(x_low) & is.finite(w_low) & w_low > 0
  keep_high <- is.finite(x_high) & is.finite(w_high) & w_high > 0
  x_low <- as.numeric(x_low[keep_low])
  w_low <- as.numeric(w_low[keep_low])
  x_high <- as.numeric(x_high[keep_high])
  w_high <- as.numeric(w_high[keep_high])

  estimate <- weighted_mean(x_high, w_high) - weighted_mean(x_low, w_low)
  n_low <- length(x_low)
  n_high <- length(x_high)
  draws <- replicate(reps, {
    id_low <- sample.int(n_low, n_low, replace = TRUE)
    id_high <- sample.int(n_high, n_high, replace = TRUE)
    weighted_mean(x_high[id_high], w_high[id_high]) -
      weighted_mean(x_low[id_low], w_low[id_low])
  })
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble(delta = estimate, delta_ci_low = ci[1], delta_ci_high = ci[2])
}

weighted_difference_bootstrap_pps <- function(
    x_reference, w_reference, x_comparison, w_comparison,
    reps = N_BOOT_NOCTURNALITY) {
  keep_reference <- is.finite(x_reference) &
    is.finite(w_reference) & w_reference > 0
  keep_comparison <- is.finite(x_comparison) &
    is.finite(w_comparison) & w_comparison > 0
  x_reference <- as.numeric(x_reference[keep_reference])
  w_reference <- as.numeric(w_reference[keep_reference])
  x_comparison <- as.numeric(x_comparison[keep_comparison])
  w_comparison <- as.numeric(w_comparison[keep_comparison])
  if (!length(x_reference) || !length(x_comparison)) {
    stop("Both groups require at least one finite, positively weighted event.")
  }

  estimate <- weighted_mean(x_comparison, w_comparison) -
    weighted_mean(x_reference, w_reference)
  n_reference <- length(x_reference)
  n_comparison <- length(x_comparison)
  draws <- replicate(reps, {
    # This implements the final Methods literally: detection records are
    # selected with probability proportional to group size, and group size
    # remains the weight in the mean calculated from each bootstrap sample.
    id_reference <- sample.int(
      n_reference,
      n_reference,
      replace = TRUE,
      prob = w_reference
    )
    id_comparison <- sample.int(
      n_comparison,
      n_comparison,
      replace = TRUE,
      prob = w_comparison
    )
    weighted_mean(
      x_comparison[id_comparison],
      w_comparison[id_comparison]
    ) - weighted_mean(
      x_reference[id_reference],
      w_reference[id_reference]
    )
  })
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble(delta = estimate, delta_ci_low = ci[1], delta_ci_high = ci[2])
}

select_extreme_cameras <- function(site_covs, covariate, n_each = N_Q_EXTREME_CAMERAS) {
  ranked <- site_covs |>
    transmute(
      camera_name = as.character(Camera.ID),
      value = suppressWarnings(as.numeric(.data[[covariate]]))
    ) |>
    filter(is.finite(value)) |>
    arrange(value, as.numeric(camera_name))

  if (nrow(ranked) < 2 * n_each) {
    stop("Fewer than ", 2 * n_each, " finite cameras for ", covariate, ".")
  }

  list(
    low = head(ranked$camera_name, n_each),
    high = tail(ranked$camera_name, n_each),
    assignments = bind_rows(
      head(ranked, n_each) |> mutate(quartile_group = "Q1"),
      tail(ranked, n_each) |> mutate(quartile_group = "Q4")
    )
  )
}

choose_delta_estimator <- function(n1, n2) {
  if (min(n1, n2) < DELTA1_THRESHOLD) "Dhat1" else "Dhat4"
}

standard_delta <- function(a, b, estimator = choose_delta_estimator(length(a), length(b))) {
  a <- as.numeric(a[is.finite(a)])
  b <- as.numeric(b[is.finite(b)])
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  suppressWarnings(tryCatch(
    as.numeric(overlap::overlapEst(a, b, type = estimator)),
    error = function(e) NA_real_
  ))
}

bootstrap_standard_delta <- function(
    a, b, estimator, reps = N_BOOT_INTERSPECIES_OVERLAP) {
  a <- as.numeric(a[is.finite(a)])
  b <- as.numeric(b[is.finite(b)])
  n1 <- length(a)
  n2 <- length(b)
  estimate <- standard_delta(a, b, estimator)
  if (reps >= 100L && N_CORES > 1L) {
    cores <- min(N_CORES, reps)
    chunk_sizes <- rep(reps %/% cores, cores)
    chunk_sizes[seq_len(reps %% cores)] <-
      chunk_sizes[seq_len(reps %% cores)] + 1L
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(overlap))
    parallel::clusterExport(
      cluster,
      c("a", "b", "n1", "n2", "estimator"),
      envir = environment()
    )
    parallel::clusterSetRNGStream(
      cluster,
      iseed = sample.int(.Machine$integer.max, 1)
    )
    draw_parts <- parallel::parLapply(cluster, chunk_sizes, function(n_reps) {
      replicate(n_reps, {
        aa <- sample(a, n1, replace = TRUE)
        bb <- sample(b, n2, replace = TRUE)
        suppressWarnings(tryCatch(
          as.numeric(overlap::overlapEst(aa, bb, type = estimator)),
          error = function(e) NA_real_
        ))
      })
    })
    draws <- unlist(draw_parts, use.names = FALSE)
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  } else {
    draws <- replicate(reps, {
      standard_delta(
        sample(a, n1, replace = TRUE),
        sample(b, n2, replace = TRUE),
        estimator
      )
    })
  }
  draws <- draws[is.finite(draws)]
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble(
    delta = estimate,
    bootstrap_se = sd(draws, na.rm = TRUE),
    ci_low = ci[1],
    ci_high = ci[2],
    n_boot_success = length(draws)
  )
}

permutation_standard_delta <- function(
    a, b, estimator, reps = N_RANDOMIZATIONS_DHAT1) {
  a <- as.numeric(a[is.finite(a)])
  b <- as.numeric(b[is.finite(b)])
  n1 <- length(a)
  pooled <- c(a, b)
  estimate <- standard_delta(a, b, estimator)

  if (reps >= 100L && N_CORES > 1L) {
    cores <- min(N_CORES, reps)
    chunk_sizes <- rep(reps %/% cores, cores)
    chunk_sizes[seq_len(reps %% cores)] <-
      chunk_sizes[seq_len(reps %% cores)] + 1L
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(overlap))
    parallel::clusterExport(
      cluster,
      c("pooled", "n1", "estimator"),
      envir = environment()
    )
    parallel::clusterSetRNGStream(
      cluster,
      iseed = sample.int(.Machine$integer.max, 1)
    )
    null_parts <- parallel::parLapply(cluster, chunk_sizes, function(n_reps) {
      replicate(n_reps, {
        id <- sample.int(length(pooled), n1, replace = FALSE)
        suppressWarnings(tryCatch(
          as.numeric(overlap::overlapEst(
            pooled[id],
            pooled[-id],
            type = estimator
          )),
          error = function(e) NA_real_
        ))
      })
    })
    null <- unlist(null_parts, use.names = FALSE)
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  } else {
    null <- replicate(reps, {
      id <- sample.int(length(pooled), n1, replace = FALSE)
      standard_delta(pooled[id], pooled[-id], estimator)
    })
  }
  null <- null[is.finite(null)]

  tibble(
    delta = estimate,
    null_mean = mean(null),
    null_se = sd(null),
    p_value = (1 + sum(null <= estimate)) / (1 + length(null)),
    n_randomization_success = length(null)
  )
}

fit_weighted_activity <- function(times, weights, reps = 1L) {
  keep <- is.finite(times) & is.finite(weights) & weights > 0
  times <- as.numeric(times[keep])
  weights <- as.numeric(weights[keep])
  if (length(times) < 5) return(NULL)
  suppressMessages(tryCatch(
    activity::fitact(
      times,
      wt = weights,
      reps = reps,
      sample = "none",
      adj = 1,
      show = FALSE
    ),
    error = function(e) NULL
  ))
}

.fast_ckern_loaded <- FALSE
.fast_ckern_max_validation_error <- NA_real_

load_fast_ckern <- function() {
  if (.fast_ckern_loaded) return(TRUE)
  source_file <- file.path(PROJECT_DIR, "src", "fast_ckern.cpp")
  if (!requireNamespace("Rcpp", quietly = TRUE) || !file.exists(source_file)) {
    return(FALSE)
  }
  loaded <- tryCatch(
    {
      suppressWarnings(
        Rcpp::sourceCpp(
          source_file,
          env = globalenv(),
          rebuild = FALSE,
          showOutput = FALSE,
          verbose = FALSE
        )
      )
      TRUE
    },
    error = function(e) {
      warning("Compiled compareCkern backend unavailable: ", conditionMessage(e))
      FALSE
    }
  )
  .fast_ckern_loaded <<- loaded
  loaded
}

compare_ckern_compiled <- function(fit1, fit2, reps) {
  if (!load_fast_ckern()) return(NULL)

  grid <- fit1@pdf[, 1]
  if (diff(range(grid)) != 2 * pi) return(NULL)
  pooled_y <- c(fit1@data, fit2@data)
  pooled_w <- c(fit1@wt, fit2@wt)
  if (length(fit1@wt) == 1L) {
    pooled_w <- c(rep(1, length(fit1@data)), fit2@wt)
  }
  if (length(fit2@wt) == 1L) {
    pooled_w <- c(
      pooled_w[seq_along(fit1@data)],
      rep(1, length(fit2@data))
    )
  }
  n1 <- length(fit1@data)
  n_total <- length(pooled_y)
  grid_step <- grid[2] - grid[1]
  left_index <- floor((pooled_y - grid[1]) / grid_step)
  left_index <- pmax(0L, pmin(as.integer(left_index), length(grid) - 2L))
  interpolation_weight <- (
    pooled_y - grid[left_index + 1L]
  ) / grid_step

  kernel <- build_vm_kernel_cpp(
    pooled_y,
    grid,
    fit1@bw * fit1@adj,
    N_CORES
  )
  batch_size <- max(1L, min(64L, floor(2e6 / n_total)))
  null <- numeric(reps)
  cursor <- 1L
  validated <- is.finite(.fast_ckern_max_validation_error)

  while (cursor <= reps) {
    n_batch <- min(batch_size, reps - cursor + 1L)
    sampled <- matrix(
      sample.int(
        n_total,
        n_batch * n_total,
        replace = TRUE,
        prob = pooled_w
      ),
      nrow = n_batch
    )
    batch_result <- ckern_overlap_batch_cpp(
      sampled,
      kernel,
      left_index,
      interpolation_weight,
      n1,
      N_CORES
    )

    if (!validated) {
      sample1 <- sampled[1, seq_len(n1)]
      sample2 <- sampled[1, (n1 + 1L):n_total]
      model1 <- activity::fitact(
        pooled_y[sample1],
        sample = "none",
        bw = fit1@bw,
        adj = fit1@adj,
        show = FALSE
      )
      model2 <- activity::fitact(
        pooled_y[sample2],
        sample = "none",
        bw = fit1@bw,
        adj = fit1@adj,
        show = FALSE
      )
      reference <- activity:::ovl4(model1, model2)
      error <- abs(reference - batch_result[1])
      .fast_ckern_max_validation_error <<- error
      if (!is.finite(error) || error > 1e-10) {
        stop(
          "Compiled compareCkern backend failed equivalence check: error = ",
          format(error, scientific = TRUE)
        )
      }
      validated <- TRUE
    }

    null[cursor:(cursor + n_batch - 1L)] <- batch_result
    cursor <- cursor + n_batch
  }

  observed <- activity:::ovl4(fit1, fit2)
  c(
    obs = observed,
    null = mean(null),
    seNull = sd(null),
    pNull = mean(null <= observed)
  )
}

compare_ckern_parallel <- function(
    fit1, fit2, reps = N_RANDOMIZATIONS_PART1) {
  compiled <- compare_ckern_compiled(fit1, fit2, reps)
  if (!is.null(compiled)) return(compiled)

  if (reps < 100L || N_CORES <= 1L) {
    return(activity::compareCkern(fit1, fit2, reps = reps))
  }

  bounds <- range(fit1@pdf[, 1])
  if (!all(bounds == range(fit2@pdf[, 1]))) {
    stop("Activity-model bounds differ.")
  }
  if (diff(bounds) == 2 * pi) bounds <- NULL

  observed <- activity:::ovl4(fit1, fit2)
  y1 <- fit1@data
  y2 <- fit2@data
  w1 <- fit1@wt
  w2 <- fit2@wt
  if (length(w1) == 1L) w1 <- rep(1, length(y1))
  if (length(w2) == 1L) w2 <- rep(1, length(y2))
  pooled_y <- c(y1, y2)
  pooled_w <- c(w1, w2)
  n1 <- length(y1)
  n_total <- length(pooled_y)
  bandwidth <- fit1@bw
  adjustment <- fit1@adj

  cores <- min(N_CORES, reps)
  chunk_sizes <- rep(reps %/% cores, cores)
  chunk_sizes[seq_len(reps %% cores)] <-
    chunk_sizes[seq_len(reps %% cores)] + 1L
  cluster <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, library(activity))
  parallel::clusterExport(
    cluster,
    c(
      "pooled_y", "pooled_w", "n1", "n_total", "bandwidth",
      "adjustment", "bounds"
    ),
    envir = environment()
  )
  parallel::clusterSetRNGStream(
    cluster,
    iseed = sample.int(.Machine$integer.max, 1)
  )
  result_parts <- parallel::parLapply(cluster, chunk_sizes, function(n_reps) {
    replicate(n_reps, {
      sampled <- sample.int(
        n_total,
        n_total,
        replace = TRUE,
        prob = pooled_w
      )
      model1 <- activity::fitact(
        pooled_y[sampled[seq_len(n1)]],
        sample = "none",
        bw = bandwidth,
        adj = adjustment,
        bounds = bounds,
        show = FALSE
      )
      model2 <- activity::fitact(
        pooled_y[sampled[(n1 + 1L):n_total]],
        sample = "none",
        bw = bandwidth,
        adj = adjustment,
        bounds = bounds,
        show = FALSE
      )
      activity:::ovl4(model1, model2)
    })
  })
  null <- unlist(result_parts, use.names = FALSE)
  parallel::stopCluster(cluster)
  on.exit(NULL, add = FALSE)
  c(
    obs = observed,
    null = mean(null),
    seNull = sd(null),
    pNull = mean(null <= observed)
  )
}

activity_curve <- function(df, weighted = TRUE, grid_n = N_ACTIVITY_GRID) {
  weights <- if (weighted) df$sighting_quantity else rep(1, nrow(df))
  fit <- fit_weighted_activity(df$time_rad, weights, reps = 1L)
  if (is.null(fit)) return(NULL)

  raw <- as.data.frame(fit@pdf)
  names(raw)[1:2] <- c("radian", "density_radian")
  hour <- (raw$radian %% (2 * pi)) * 24 / (2 * pi)
  density_hour <- raw$density_radian * (2 * pi / 24)
  base <- tibble(hour = hour, density = density_hour) |>
    arrange(hour) |>
    group_by(hour) |>
    summarise(density = mean(density), .groups = "drop")

  grid <- seq(0, 24, length.out = grid_n)
  periodic <- bind_rows(
    base |> mutate(hour = hour - 24),
    base,
    base |> mutate(hour = hour + 24)
  ) |>
    arrange(hour)
  y <- approx(periodic$hour, periodic$density, xout = grid, rule = 2)$y
  y <- normalize_curve(grid, pmax(y, 0))
  tibble(hour = grid, density = y)
}

curve_overlap <- function(curve1, curve2) {
  if (is.null(curve1) || is.null(curve2)) return(NA_real_)
  stopifnot(identical(curve1$hour, curve2$hour))
  trapz(curve1$hour, pmin(curve1$density, curve2$density))
}

adaptive_overlap_test <- function(
    a, b, w_a = NULL, w_b = NULL,
    dhat4_reps = N_RANDOMIZATIONS_PART1,
    dhat1_reps = N_RANDOMIZATIONS_DHAT1,
    bootstrap_reps = N_BOOT_DHAT1) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (is.null(w_a)) w_a <- rep(1, length(a))
  if (is.null(w_b)) w_b <- rep(1, length(b))
  keep_a <- is.finite(a) & is.finite(w_a) & w_a > 0
  keep_b <- is.finite(b) & is.finite(w_b) & w_b > 0
  a <- a[keep_a]
  b <- b[keep_b]
  w_a <- as.numeric(w_a[keep_a])
  w_b <- as.numeric(w_b[keep_b])
  estimator <- choose_delta_estimator(length(a), length(b))

  if (estimator == "Dhat4") {
    fit_a <- fit_weighted_activity(a, w_a, reps = 1L)
    fit_b <- fit_weighted_activity(b, w_b, reps = 1L)
    if (is.null(fit_a) || is.null(fit_b)) stop("Weighted activity model failed.")
    ck <- compare_ckern_parallel(fit_a, fit_b, reps = dhat4_reps)
    return(tibble(
      n1 = length(a),
      n2 = length(b),
      n_individuals1 = sum(w_a),
      n_individuals2 = sum(w_b),
      estimator = estimator,
      delta = unname(as.numeric(ck["obs"])),
      bootstrap_se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      null_mean = unname(as.numeric(ck["null"])),
      null_se = unname(as.numeric(ck["seNull"])),
      p_value = unname(as.numeric(ck["pNull"])),
      test = "weighted compareCkern",
      n_randomizations = dhat4_reps
    ))
  }

  # Dhat1 is required only for low-event groups. overlapEst() has no weight
  # argument; the detection event remains the independent sampling unit.
  perm <- permutation_standard_delta(a, b, estimator, dhat1_reps)
  boot <- bootstrap_standard_delta(a, b, estimator, bootstrap_reps)
  tibble(
    n1 = length(a),
    n2 = length(b),
    n_individuals1 = sum(w_a),
    n_individuals2 = sum(w_b),
    estimator = estimator,
    delta = perm$delta,
    bootstrap_se = boot$bootstrap_se,
    ci_low = boot$ci_low,
    ci_high = boot$ci_high,
    null_mean = perm$null_mean,
    null_se = perm$null_se,
    p_value = perm$p_value,
    test = "event-level Dhat1 permutation",
    n_randomizations = perm$n_randomization_success
  )
}

daypart_layers <- function() {
  list(
    annotate("rect", xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf,
             fill = "grey55", alpha = 0.20),
    annotate("rect", xmin = 5, xmax = 7, ymin = -Inf, ymax = Inf,
             fill = "#CFE8FF", alpha = 0.45),
    annotate("rect", xmin = 7, xmax = 17, ymin = -Inf, ymax = Inf,
             fill = "#FFF4B3", alpha = 0.38),
    annotate("rect", xmin = 17, xmax = 19, ymin = -Inf, ymax = Inf,
             fill = "#CFE8FF", alpha = 0.45),
    annotate("rect", xmin = 19, xmax = 24, ymin = -Inf, ymax = Inf,
             fill = "grey55", alpha = 0.20)
  )
}

save_figure <- function(plot, path_no_extension, width, height, dpi = 300) {
  ggsave(
    paste0(path_no_extension, ".png"),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    paste0(path_no_extension, ".pdf"),
    plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

format_p <- function(p) {
  case_when(
    !is.finite(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}
