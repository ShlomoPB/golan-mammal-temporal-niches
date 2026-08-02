load_analysis_data <- function() {
  assert_file(DATA_CAPTURES)
  assert_file(DATA_SITECOVS)

  captures_md5 <- tolower(unname(tools::md5sum(DATA_CAPTURES)))
  if (!identical(captures_md5, EXPECTED_CAPTURES_MD5)) {
    stop(
      "Capture input checksum mismatch. Expected the camera-local-time file ",
      EXPECTED_CAPTURES_MD5, " but found ", captures_md5,
      ". Do not use a file that converts Camelot timestamps from UTC."
    )
  }

  site_covs <- read.csv(DATA_SITECOVS, check.names = FALSE) |>
    as_tibble() |>
    mutate(
      camera_name = as.character(Camera.ID),
      land_use = recode(
        Polygon.Type,
        "High Culling" = "High culling",
        "Low Culling" = "Low culling",
        "Nature Reserve" = "Nature reserves",
        "Minefields" = "Minefields"
      ),
      land_use = factor(land_use, levels = LAND_USE_LEVELS)
    )

  if (nrow(site_covs) != 60L || n_distinct(site_covs$camera_name) != 60L) {
    stop("SiteCovs must contain exactly 60 unique camera stations.")
  }
  if (!all(table(site_covs$land_use) == 15L)) {
    stop("Expected exactly 15 cameras in each land-use category.")
  }

  captures <- read.csv(DATA_CAPTURES, check.names = FALSE) |>
    clean_names() |>
    as_tibble() |>
    mutate(
      camera_name = as.character(camera_name),
      species = trimws(as.character(species)),
      time_rad = suppressWarnings(as.numeric(time_rad)),
      sighting_quantity = suppressWarnings(as.numeric(sighting_quantity)),
      sighting_quantity = if_else(
        is.finite(sighting_quantity) & sighting_quantity > 0,
        sighting_quantity,
        1
      )
    ) |>
    left_join(
      site_covs |> select(camera_name, land_use),
      by = "camera_name"
    ) |>
    mutate(
      solar_hour = (time_rad %% (2 * pi)) * 24 / (2 * pi),
      nocturnality_h = pmin(abs(solar_hour - 12), 24 - abs(solar_hour - 12))
    )

  timestamp_sentinel <- captures |>
    filter(
      camera_name == "1",
      trap_station_session_id == 7,
      species == "Bird unknown",
      date_time == "2020-09-11 17:09:05"
    )
  if (nrow(timestamp_sentinel) != 1L) {
    stop(
      "Camera-local timestamp validation failed. The canonical Camelot event ",
      "must remain 2020-09-11 17:09:05, not be shifted to 20:09:05."
    )
  }

  focal_raw <- captures |>
    filter(
      species %in% SPECIES,
      is.finite(time_rad),
      is.finite(sighting_quantity),
      sighting_quantity > 0,
      !is.na(land_use)
    ) |>
    mutate(
      species = factor(species, levels = SPECIES),
      land_use = factor(land_use, levels = LAND_USE_LEVELS)
    )

  # Camelot can repeat an event when adjacent deployment sessions meet. Exact
  # camera-species-timestamp repeats are one independent detection under the
  # manuscript's 30-minute rule and must be collapsed consistently in every
  # analysis, not only when reporting Q1/Q4 sample sizes.
  focal_exact_duplicate_rows <- focal_raw |>
    group_by(camera_name, species, date_time) |>
    filter(n() > 1L, n_distinct(trap_station_session_id) > 1L) |>
    mutate(export_rows_for_event = n()) |>
    ungroup() |>
    arrange(camera_name, species, date_time, trap_station_session_id)

  focal <- focal_raw |>
    arrange(camera_name, species, date_time, trap_station_session_id) |>
    distinct(camera_name, species, date_time, .keep_all = TRUE)

  list(
    captures = captures,
    focal = focal,
    site_covs = site_covs,
    focal_exact_duplicate_rows = focal_exact_duplicate_rows
  )
}
