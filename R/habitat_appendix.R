# Appendix C habitat-confounding analysis.
#
# Camera-centered land-cover proportions are descriptive at the camera level.
# Management-category tests use the 12 sampling polygons as independent units
# after averaging the five cameras in each polygon. Exact P values enumerate
# every allocation of 12 polygons to four labelled groups of three.

habitat_management_levels <- c(
  "High-culling pasture",
  "Low-culling pasture",
  "Nature reserve",
  "Minefield"
)

build_balanced_habitat_assignments <- function(
  n_units = 12L,
  group_size = 3L
) {
  triples <- combn(n_units, group_size, simplify = FALSE)
  triple_key <- vapply(
    triples,
    function(x) paste(sort(x), collapse = "-"),
    character(1)
  )
  triple_lookup <- setNames(seq_along(triples), triple_key)
  n_assignments <- choose(12, 3) * choose(9, 3) * choose(6, 3)
  labels <- matrix(4L, nrow = n_assignments, ncol = n_units)
  combo_ids <- matrix(NA_integer_, nrow = n_assignments, ncol = 4L)
  row_index <- 0L

  for (g1 in triples) {
    remaining_1 <- setdiff(seq_len(n_units), g1)
    for (g2 in combn(remaining_1, group_size, simplify = FALSE)) {
      remaining_2 <- setdiff(remaining_1, g2)
      for (g3 in combn(remaining_2, group_size, simplify = FALSE)) {
        g4 <- setdiff(remaining_2, g3)
        row_index <- row_index + 1L
        labels[row_index, g1] <- 1L
        labels[row_index, g2] <- 2L
        labels[row_index, g3] <- 3L
        combo_ids[row_index, ] <- c(
          triple_lookup[[paste(sort(g1), collapse = "-")]],
          triple_lookup[[paste(sort(g2), collapse = "-")]],
          triple_lookup[[paste(sort(g3), collapse = "-")]],
          triple_lookup[[paste(sort(g4), collapse = "-")]]
        )
      }
    }
  }

  if (row_index != n_assignments) {
    stop("Exact habitat assignment enumeration is incomplete.")
  }
  list(
    labels = labels,
    combo_ids = combo_ids,
    triples = triples,
    triple_lookup = triple_lookup,
    n_assignments = n_assignments
  )
}

exact_habitat_test <- function(response, group, assignments) {
  group <- factor(group, levels = habitat_management_levels)
  group_codes <- as.integer(group)
  if (length(response) != 12L || !all(table(group_codes) == 3L)) {
    stop("Habitat tests require 12 polygons and three polygons per category.")
  }

  distance_squared <- as.matrix(stats::dist(response))^2
  total_ss <- sum(distance_squared[upper.tri(distance_squared)]) /
    length(response)
  triple_wss <- vapply(assignments$triples, function(index) {
    within <- distance_squared[index, index, drop = FALSE]
    sum(within[upper.tri(within)]) / length(index)
  }, numeric(1))
  permutation_within_ss <- rowSums(matrix(
    triple_wss[assignments$combo_ids],
    nrow = assignments$n_assignments,
    ncol = 4L
  ))
  permutation_between_ss <- total_ss - permutation_within_ss
  permutation_f <- (permutation_between_ss / 3) /
    (permutation_within_ss / 8)

  observed_combo_ids <- vapply(seq_len(4L), function(code) {
    key <- paste(which(group_codes == code), collapse = "-")
    assignments$triple_lookup[[key]]
  }, integer(1))
  observed_within_ss <- sum(triple_wss[observed_combo_ids])
  observed_between_ss <- total_ss - observed_within_ss
  observed_f <- (observed_between_ss / 3) / (observed_within_ss / 8)

  list(
    f = observed_f,
    eta_squared = observed_between_ss / total_ss,
    p = mean(permutation_f >= observed_f - 1e-12),
    n_assignments = assignments$n_assignments
  )
}

extract_habitat_buffer_composition <- function(
  cameras,
  landcover,
  valid_codes,
  radius_m
) {
  buffers <- sf::st_buffer(cameras, dist = radius_m, nQuadSegs = 90)
  buffer_area_m2 <- as.numeric(sf::st_area(buffers))
  extracted <- terra::extract(
    landcover,
    terra::vect(buffers),
    cells = TRUE,
    exact = TRUE
  ) |>
    tibble::as_tibble()
  names(extracted)[2] <- "code"
  if (!"fraction" %in% names(extracted)) {
    stop("Exact habitat extraction did not return raster-cell fractions.")
  }

  cell_area_m2 <- prod(terra::res(landcover))
  extracted <- extracted |>
    dplyr::mutate(
      code = as.integer(code),
      intersected_area_m2 = fraction * cell_area_m2,
      valid_landcover = !is.na(code) & code != -128L
    )

  coverage <- extracted |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      valid_raster_area_m2 = sum(
        intersected_area_m2[valid_landcover],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      buffer_area_m2 = buffer_area_m2[ID],
      valid_raster_coverage_pct =
        100 * valid_raster_area_m2 / buffer_area_m2
    )
  if (any(coverage$valid_raster_area_m2 <= 0)) {
    stop("At least one camera buffer has no valid land-cover cells.")
  }

  class_wide <- extracted |>
    dplyr::filter(valid_landcover) |>
    dplyr::group_by(ID, code) |>
    dplyr::summarise(
      area_m2 = sum(intersected_area_m2),
      .groups = "drop"
    ) |>
    tidyr::complete(
      ID = seq_len(nrow(cameras)),
      code = valid_codes,
      fill = list(area_m2 = 0)
    ) |>
    dplyr::left_join(
      coverage |>
        dplyr::select(ID, valid_raster_area_m2),
      by = "ID"
    ) |>
    dplyr::mutate(
      proportion = area_m2 / valid_raster_area_m2,
      column = paste0("prop_class_", code)
    ) |>
    dplyr::select(ID, column, proportion) |>
    tidyr::pivot_wider(
      names_from = column,
      values_from = proportion,
      values_fill = 0
    )

  result <- sf::st_drop_geometry(cameras) |>
    dplyr::mutate(
      ID = dplyr::row_number(),
      buffer_radius_m = as.integer(radius_m)
    ) |>
    dplyr::left_join(coverage, by = "ID") |>
    dplyr::left_join(class_wide, by = "ID") |>
    dplyr::mutate(
      natural_open = prop_class_31 + prop_class_32,
      natural_woody = prop_class_33 + prop_class_34,
      structural_woody = natural_woody + prop_class_35
    ) |>
    dplyr::select(-ID)

  class_columns <- grep("^prop_class_", names(result), value = TRUE)
  if (any(abs(rowSums(result[, class_columns]) - 1) > 1e-8)) {
    stop("Land-cover class proportions do not sum to one.")
  }
  result
}

run_habitat_appendix <- function() {
  message("Running Appendix C habitat-confounding analysis")
  required_files <- c(
    HABITAT_RASTER,
    HABITAT_CAMERAS,
    HABITAT_POLYGONS,
    HABITAT_KEY,
    HABITAT_FIELD_VEGETATION
  )
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop(
      "Missing Appendix C habitat input(s): ",
      paste(missing_files, collapse = ", ")
    )
  }

  cameras <- sf::st_read(HABITAT_CAMERAS, quiet = TRUE)
  polygons <- sf::st_read(HABITAT_POLYGONS, quiet = TRUE) |>
    sf::st_make_valid()
  key <- readr::read_csv(HABITAT_KEY, show_col_types = FALSE)
  field <- readr::read_csv(
    HABITAT_FIELD_VEGETATION,
    show_col_types = FALSE
  )
  landcover <- terra::rast(HABITAT_RASTER)
  levels(landcover) <- NULL

  if (
    nrow(cameras) != 60L ||
      dplyr::n_distinct(cameras$camera_id) != 60L ||
      !setequal(cameras$camera_id, 1:60)
  ) {
    stop("Appendix C requires exactly 60 unique camera IDs 1-60.")
  }
  if (
    nrow(polygons) != 12L ||
      dplyr::n_distinct(polygons$polygon_id) != 12L
  ) {
    stop("Appendix C requires exactly 12 unique sampling polygons.")
  }
  if (!identical(sf::st_crs(cameras), sf::st_crs(polygons))) {
    stop("Camera and polygon coordinate reference systems differ.")
  }
  if (!isTRUE(sf::st_crs(cameras) == sf::st_crs(terra::crs(landcover)))) {
    stop("Vector and land-cover raster coordinate reference systems differ.")
  }
  if (!all(terra::res(landcover) == c(25, 25))) {
    stop("The Appendix C land-cover raster must retain its 25 m cells.")
  }

  spatial_membership <- sf::st_join(
    cameras |>
      dplyr::select(camera_id),
    polygons |>
      dplyr::select(
        spatial_polygon_id = polygon_id,
        spatial_management = management_category
      ),
    join = sf::st_within,
    left = TRUE
  ) |>
    sf::st_drop_geometry()
  if (
    nrow(spatial_membership) != 60L ||
      any(is.na(spatial_membership$spatial_polygon_id))
  ) {
    stop("At least one camera does not fall in exactly one sampling polygon.")
  }
  membership_check <- sf::st_drop_geometry(cameras) |>
    dplyr::left_join(spatial_membership, by = "camera_id")
  if (
    any(membership_check$polygon_id != membership_check$spatial_polygon_id) ||
      any(
        membership_check$management_category !=
          membership_check$spatial_management
      )
  ) {
    stop("Camera attributes disagree with spatial polygon membership.")
  }

  design <- sf::st_drop_geometry(cameras) |>
    dplyr::count(
      management_category,
      polygon_id,
      sampling_polygon,
      name = "n_cameras"
    )
  if (
    nrow(design) != 12L ||
      any(design$n_cameras != 5L) ||
      any(table(design$management_category) != 3L)
  ) {
    stop("Expected five cameras per polygon and three polygons per category.")
  }

  expected_codes <- c(-128L, 11L, 12L, 21L, 22L, 23L, 24L, 25L,
                      31L, 32L, 33L, 34L, 35L, 41L, 42L, 51L, 52L)
  if (!setequal(key$code, expected_codes)) {
    stop("The land-cover lookup does not contain the expected 16 classes and NoData.")
  }
  valid_codes <- sort(key$code[key$valid_landcover])

  radii <- c(250L, 500L, 1000L)
  camera_by_scale <- lapply(radii, function(radius) {
    extract_habitat_buffer_composition(
      cameras,
      landcover,
      valid_codes,
      radius
    )
  })
  names(camera_by_scale) <- as.character(radii)
  camera_all <- dplyr::bind_rows(camera_by_scale)

  for (radius in radii) {
    readr::write_csv(
      camera_by_scale[[as.character(radius)]],
      file.path(
        OUTPUT_DIRS$diagnostics,
        paste0("habitat_camera_landcover_", radius, "m.csv")
      )
    )
  }

  class_columns <- grep("^prop_class_", names(camera_all), value = TRUE)
  polygon_means <- camera_all |>
    dplyr::group_by(
      buffer_radius_m,
      management_category,
      polygon_id,
      sampling_polygon
    ) |>
    dplyr::summarise(
      n_cameras = dplyr::n(),
      mean_valid_raster_coverage_pct = mean(valid_raster_coverage_pct),
      dplyr::across(
        dplyr::all_of(c(
          class_columns,
          "natural_open",
          "natural_woody",
          "structural_woody"
        )),
        mean
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      buffer_radius_m,
      factor(management_category, levels = habitat_management_levels),
      polygon_id
    )
  if (any(polygon_means$n_cameras != 5L)) {
    stop("Not every polygon mean contains five camera values.")
  }
  readr::write_csv(
    polygon_means,
    file.path(OUTPUT_DIRS$diagnostics, "habitat_polygon_means.csv")
  )

  assignments <- build_balanced_habitat_assignments()
  tests <- dplyr::bind_rows(lapply(radii, function(radius) {
    scale_data <- polygon_means |>
      dplyr::filter(buffer_radius_m == radius) |>
      dplyr::arrange(polygon_id)
    dplyr::bind_rows(lapply(
      c("structural_woody", "natural_open"),
      function(variable) {
        test <- exact_habitat_test(
          scale_data[[variable]],
          scale_data$management_category,
          assignments
        )
        tibble::tibble(
          buffer_radius_m = radius,
          variable = variable,
          f_statistic = test$f,
          eta_squared = test$eta_squared,
          p_exact = test$p,
          exact_balanced_assignments = test$n_assignments
        )
      }
    ))
  }))
  readr::write_csv(
    tests,
    file.path(OUTPUT_DIRS$additional_results, "habitat_exact_tests.csv")
  )

  primary_camera <- camera_by_scale[["500"]]
  table_c1 <- dplyr::bind_rows(lapply(
    c("structural_woody", "natural_open"),
    function(variable) {
      primary_camera |>
        dplyr::group_by(management_category) |>
        dplyr::summarise(
          variable = variable,
          n_cameras = dplyr::n(),
          mean_proportion = mean(.data[[variable]]),
          sd_proportion = stats::sd(.data[[variable]]),
          median_proportion = stats::median(.data[[variable]]),
          minimum_proportion = min(.data[[variable]]),
          maximum_proportion = max(.data[[variable]]),
          .groups = "drop"
        )
    }
  )) |>
    dplyr::mutate(
      management_category = factor(
        management_category,
        levels = habitat_management_levels
      ),
      variable_label = dplyr::recode(
        variable,
        structural_woody = "Structural woody cover",
        natural_open = "Natural open/herbaceous cover"
      ),
      dplyr::across(
        dplyr::ends_with("_proportion"),
        ~100 * .x,
        .names = "{sub('_proportion$', '_percent', .col)}"
      )
    ) |>
    dplyr::arrange(variable, management_category) |>
    dplyr::select(
      management_category,
      variable,
      variable_label,
      n_cameras,
      dplyr::ends_with("_proportion"),
      dplyr::ends_with("_percent")
    )
  readr::write_csv(
    table_c1,
    file.path(
      OUTPUT_DIRS$supplementary_tables,
      "Table_C1_habitat_composition_500m.csv"
    )
  )

  field_per_camera <- field |>
    dplyr::distinct(camera_id, vegetation_structure)
  if (
    nrow(field_per_camera) != 60L ||
      any(dplyr::count(field_per_camera, camera_id)$n != 1L)
  ) {
    stop("Each camera must have one field vegetation classification.")
  }
  field_joined <- field_per_camera |>
    dplyr::left_join(
      sf::st_drop_geometry(cameras) |>
        dplyr::select(
          camera_id,
          polygon_id,
          sampling_polygon,
          management_category
        ),
      by = "camera_id"
    )
  field_summary <- field_joined |>
    dplyr::count(
      management_category,
      vegetation_structure,
      name = "n_cameras"
    ) |>
    tidyr::complete(
      management_category = habitat_management_levels,
      vegetation_structure = sort(unique(field_joined$vegetation_structure)),
      fill = list(n_cameras = 0)
    ) |>
    dplyr::group_by(management_category) |>
    dplyr::mutate(
      percent_within_management = 100 * n_cameras / sum(n_cameras)
    ) |>
    dplyr::ungroup()
  readr::write_csv(
    field_summary,
    file.path(
      OUTPUT_DIRS$diagnostics,
      "field_vegetation_structure_summary.csv"
    )
  )

  polygon_order <- sf::st_drop_geometry(polygons) |>
    dplyr::arrange(polygon_id) |>
    dplyr::select(polygon_id, management_category)
  vegetation_levels <- sort(unique(field_joined$vegetation_structure))
  field_matrix <- field_joined |>
    dplyr::count(polygon_id, vegetation_structure, name = "n") |>
    tidyr::complete(
      polygon_id = polygon_order$polygon_id,
      vegetation_structure = vegetation_levels,
      fill = list(n = 0)
    ) |>
    dplyr::arrange(
      polygon_id,
      match(vegetation_structure, vegetation_levels)
    ) |>
    tidyr::pivot_wider(names_from = vegetation_structure, values_from = n) |>
    dplyr::arrange(polygon_id) |>
    dplyr::select(dplyr::all_of(vegetation_levels)) |>
    as.matrix()
  observed_codes <- as.integer(factor(
    polygon_order$management_category,
    levels = habitat_management_levels
  ))
  chisq_statistic <- function(labels) {
    group_counts <- rowsum(field_matrix, labels, reorder = FALSE)
    expected <- outer(rowSums(group_counts), colSums(group_counts)) /
      sum(group_counts)
    sum((group_counts - expected)^2 / expected)
  }
  field_chisq <- chisq_statistic(observed_codes)
  expected_per_group <- colSums(field_matrix) / 4
  permuted_chisq <- numeric(assignments$n_assignments)
  for (group_code in seq_len(4L)) {
    permuted_counts <-
      (assignments$labels == group_code) %*% field_matrix
    permuted_chisq <- permuted_chisq + rowSums(
      sweep(permuted_counts, 2, expected_per_group, "-")^2 /
        matrix(
          expected_per_group,
          nrow = nrow(permuted_counts),
          ncol = ncol(permuted_counts),
          byrow = TRUE
        )
    )
  }
  field_test <- tibble::tibble(
    chi_square = field_chisq,
    cramers_v = sqrt(
      field_chisq /
        (60 * min(3, length(vegetation_levels) - 1))
    ),
    p_exact = mean(permuted_chisq >= field_chisq - 1e-12),
    exact_balanced_assignments = assignments$n_assignments,
    independent_polygons = 12L
  )
  readr::write_csv(
    field_test,
    file.path(
      OUTPUT_DIRS$additional_results,
      "field_vegetation_structure_exact_test.csv"
    )
  )

  primary_means <- table_c1 |>
    dplyr::select(management_category, variable, mean_proportion) |>
    tidyr::pivot_wider(names_from = management_category, values_from = mean_proportion)
  expected_means <- tibble::tribble(
    ~variable, ~`High-culling pasture`, ~`Low-culling pasture`, ~`Nature reserve`, ~Minefield,
    "structural_woody", 0.05457614006115687, 0.07964942699931819, 0.6620832408916846, 0.12583592128733015,
    "natural_open", 0.8999053798008245, 0.8226060299177546, 0.28790074123815046, 0.8116049766836773
  )
  comparison <- dplyr::left_join(
    primary_means,
    expected_means,
    by = "variable",
    suffix = c("_observed", "_expected")
  )
  for (management in habitat_management_levels) {
    difference <- abs(
      comparison[[paste0(management, "_observed")]] -
        comparison[[paste0(management, "_expected")]]
    )
    if (any(difference > 1e-10)) {
      stop("Appendix C habitat means do not match the submitted results.")
    }
  }

  list(
    n_cameras = nrow(cameras),
    n_polygons = nrow(polygons),
    camera_composition = camera_all,
    polygon_means = polygon_means,
    tests = tests,
    table_c1 = table_c1,
    field_summary = field_summary,
    field_test = field_test
  )
}
