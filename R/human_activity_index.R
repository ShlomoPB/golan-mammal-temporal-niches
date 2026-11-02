run_human_activity_index <- function(
  refit_model = FALSE,
  refit_candidates = FALSE
) {
  if (
    (isTRUE(refit_model) || isTRUE(refit_candidates)) &&
      !requireNamespace("unmarked", quietly = TRUE)
  ) {
    stop("Package 'unmarked' is required. Restore the renv environment first.")
  }

  inputs <- c(
    human_counts = DATA_HUMAN_COUNTS,
    observation_covariates = DATA_HUMAN_OBSCOVS,
    site_covariates = DATA_HUMAN_SITECOVS
  )
  expected_md5 <- c(
    human_counts = "e8ebd73e2f547d63c715310523f43628",
    observation_covariates = "c0c4e9ce9afdfdc10a51530d4292ebb2",
    site_covariates = "4829422fcd6b0637169a048616ca6869"
  )
  if (!all(file.exists(inputs))) {
    stop("One or more human-activity model inputs are missing.")
  }
  observed_md5 <- unname(tools::md5sum(inputs))
  if (!identical(observed_md5, unname(expected_md5))) {
    stop("Human-activity input checksum validation failed.")
  }

  counts <- read.csv(DATA_HUMAN_COUNTS, check.names = FALSE)
  site_covs <- read.csv(DATA_HUMAN_SITECOVS, check.names = FALSE)
  obs_covs <- read.csv(DATA_HUMAN_OBSCOVS, check.names = TRUE)

  if (!identical(dim(counts), c(60L, 7L))) {
    stop("Expected a 60-camera by 7-occasion human count history.")
  }
  if (!identical(site_covs$Camera.ID, seq_len(60L))) {
    stop("Site-covariate rows are not ordered by Camera.ID 1-60.")
  }
  expected_camera <- rep(seq_len(60L), each = 7L)
  expected_round <- rep(seq_len(7L), times = 60L)
  if (
    !identical(obs_covs$Camera.Number, expected_camera) ||
      !identical(obs_covs$Round.Number, expected_round)
  ) {
    stop("Observation covariates are not in camera-by-round order.")
  }

  site_covs$Polygon.Name <- factor(site_covs$Polygon.Name)
  site_covs$Polygon.Type <- factor(site_covs$Polygon.Type)
  obs_model <- obs_covs[
    setdiff(
      names(obs_covs),
      c("Start.date", "End.Date", "Camera.Number", "Round.Number")
    )
  ]
  factor_fields <- c(
    "Camera.Model",
    "Camera.Sensitivity",
    "Capture.Delay",
    "Vegetation.Structure",
    "Path.Type",
    "Area.Behind.Camera",
    "Bottleneck.",
    "Livestock.Present.",
    "Agricultural.Presence"
  )
  obs_model[factor_fields] <- lapply(obs_model[factor_fields], factor)
  obs_model$Days.Operational <- as.numeric(obs_model$Days.Operational)

  # Coefficients from the archived converged fit are used for the default
  # reconstruction and as starting values when the likelihood is refitted.
  final_starts <- c(
    4.70654803603723,
    -2.59421091636092,
    -3.19061655118606,
    -1.2015614192934,
    -1.43893382494163,
    1.32549208257991,
    0.879439414968176,
    0.625217931929671,
    -0.656050468552745,
    -0.365596245920933,
    2.30378090855244,
    0.865831584703568,
    -0.210893279791971,
    0.363350014370663,
    0.234819172059341
  )
  final_formula <-
    ~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1
  frame <- NULL
  final_model <- NULL
  fitted_coefficients <- final_starts
  if (isTRUE(refit_model) || isTRUE(refit_candidates)) {
    message("Refitting the final human-activity N-mixture model")
    frame <- unmarked::unmarkedFramePCount(
      y = as.matrix(counts),
      siteCovs = site_covs,
      obsCovs = obs_model
    )
    final_model <- unmarked::pcount(
      final_formula,
      frame,
      starts = final_starts
    )
    fitted_coefficients <- stats::coef(final_model)
  }

  detection_fields <- c(
    "Path.Type",
    "Camera.Model",
    "Vegetation.Structure",
    "Area.Behind.Camera"
  )
  complete_detection_covariates <- stats::complete.cases(
    obs_model[detection_fields]
  )
  detection_matrix <- stats::model.matrix(
    ~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera,
    data = obs_model[complete_detection_covariates, , drop = FALSE]
  )
  if (ncol(detection_matrix) != length(fitted_coefficients) - 1L) {
    stop("Detection design matrix does not match the fitted coefficients.")
  }
  lambda <- exp(fitted_coefficients[1])
  detection_probability_long <- rep(NA_real_, nrow(obs_model))
  detection_probability_long[complete_detection_covariates] <- stats::plogis(
    detection_matrix %*% fitted_coefficients[-1]
  )
  detection_probability <- matrix(
    detection_probability_long,
    nrow = 60L,
    ncol = 7L,
    byrow = TRUE
  )
  count_matrix <- as.matrix(counts)
  if (any(is.na(detection_probability[!is.na(count_matrix)]))) {
    stop("An observed human count is missing a required detection covariate.")
  }
  message("Calculating conditional site-level human-activity estimates")
  latent_n <- 0:(max(count_matrix, na.rm = TRUE) + 100L)
  estimate <- numeric(60L)
  ci <- matrix(NA_real_, nrow = 60L, ncol = 2L)
  for (site in seq_len(60L)) {
    log_posterior <- stats::dpois(latent_n, lambda, log = TRUE)
    observed_rounds <- which(!is.na(count_matrix[site, ]))
    for (round in observed_rounds) {
      log_posterior <- log_posterior + stats::dbinom(
        count_matrix[site, round],
        size = latent_n,
        prob = detection_probability[site, round],
        log = TRUE
      )
    }
    posterior <- exp(log_posterior - max(log_posterior))
    posterior <- posterior / sum(posterior)
    posterior_cdf <- cumsum(posterior)
    estimate[site] <- sum(latent_n * posterior)
    ci[site, 1] <- latent_n[which(posterior_cdf >= 0.025)[1]]
    ci[site, 2] <- latent_n[which(posterior_cdf >= 0.975)[1]]
  }
  scaling_max <- max(estimate, na.rm = TRUE)
  site_estimates <- data.frame(
    site = seq_along(estimate),
    Camera.ID = site_covs$Camera.ID,
    estimate = estimate,
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    human_activity_index = estimate / scaling_max,
    standardized_CI_lower = ci[, 1] / scaling_max,
    standardized_CI_upper = ci[, 2] / scaling_max
  )

  current_site_covs <- read.csv(DATA_SITECOVS, check.names = TRUE)
  if (!identical(current_site_covs$Camera.ID, site_estimates$Camera.ID)) {
    stop("Current and model-input Camera.ID order differs.")
  }
  maximum_index_difference <- max(
    abs(
      site_estimates$human_activity_index -
        current_site_covs$human.abundance
    )
  )
  if (maximum_index_difference > 5e-8) {
    stop(
      "Reconstructed human-activity index differs from SiteCovs_280825.csv; ",
      "maximum absolute difference = ",
      format(maximum_index_difference, scientific = TRUE)
    )
  }

  readr::write_csv(
    site_estimates,
    file.path(
      OUTPUT_DIRS$additional_results,
      "site_level_human_activity_estimates.csv"
    )
  )

  candidate_results <- NULL
  if (isTRUE(refit_candidates)) {
    fit_model <- function(formula) unmarked::pcount(formula, frame)
    individual_formulas <- list(
      Null = ~1 ~ 1,
      Cattle.Sightings = ~Cattle.Sightings ~ 1,
      Camera.Model = ~Camera.Model ~ 1,
      Camera.Sensitivity = ~Camera.Sensitivity ~ 1,
      Capture.Delay = ~Capture.Delay ~ 1,
      Vegetation.Structure = ~Vegetation.Structure ~ 1,
      Path.Type = ~Path.Type ~ 1,
      Area.Behind.Camera = ~Area.Behind.Camera ~ 1,
      Bottleneck = ~Bottleneck. ~ 1,
      Livestock.Present = ~Livestock.Present. ~ 1,
      Agricultural.Presence = ~Agricultural.Presence ~ 1
    )
    multivariable_formulas <- list(
      Path_Camera = ~Path.Type + Camera.Model ~ 1,
      Path_Camera_Passability =
        ~Path.Type + Camera.Model + Area.Behind.Camera ~ 1,
      Path_Camera_Passability_Bottleneck =
        ~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~ 1,
      Final =
        ~Path.Type + Camera.Model + Vegetation.Structure +
          Area.Behind.Camera ~ 1,
      Path_Camera_Vegetation =
        ~Path.Type + Camera.Model + Vegetation.Structure ~ 1
    )
    candidate_formulas <- c(individual_formulas, multivariable_formulas)
    candidate_models <- lapply(
      candidate_formulas[names(candidate_formulas) != "Final"],
      fit_model
    )
    candidate_models$Final <- final_model
    candidate_models <- candidate_models[names(candidate_formulas)]
    candidate_results <- as(
      do.call(unmarked::fitList, candidate_models),
      "data.frame"
    )
    candidate_results$model_id <- rownames(candidate_results)
    rownames(candidate_results) <- NULL
    candidate_results <- candidate_results[
      c("model_id", setdiff(names(candidate_results), "model_id"))
    ]
    readr::write_csv(
      candidate_results,
      file.path(
        OUTPUT_DIRS$additional_results,
        "human_detection_candidate_models.csv"
      )
    )
  }

  message(
    "Human-activity index reproduced; maximum absolute difference = ",
    format(maximum_index_difference, scientific = TRUE)
  )
  invisible(list(
    model = final_model,
    site_estimates = site_estimates,
    candidate_results = candidate_results,
    maximum_index_difference = maximum_index_difference
  ))
}
