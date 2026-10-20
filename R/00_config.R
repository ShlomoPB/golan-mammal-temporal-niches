options(stringsAsFactors = FALSE)

PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DATA_CAPTURES <- file.path(
  PROJECT_DIR,
  "data",
  "Independent_Captures_Sep_2025_solar.csv"
)
DATA_SITECOVS <- file.path(PROJECT_DIR, "data", "SiteCovs_280825.csv")

# Camelot Date/Time values are camera-local Asia/Jerusalem times. They must be
# assigned that timezone without conversion from UTC. This checksum guards
# against the previously used file that added the UTC+2/UTC+3 offset.
EXPECTED_CAPTURES_MD5 <- "f54e92df14a6bcb1730ecb5474370db4"
CAPTURE_TIME_BASIS <- "Camelot camera-local time (Asia/Jerusalem)"

OUTPUT_ROOT <- file.path(PROJECT_DIR, "manuscript_outputs")
OUTPUT_DIRS <- list(
  main_figures = file.path(OUTPUT_ROOT, "main_figures"),
  main_tables = file.path(OUTPUT_ROOT, "main_tables"),
  supplementary_figures = file.path(OUTPUT_ROOT, "supplementary_figures"),
  supplementary_tables = file.path(OUTPUT_ROOT, "supplementary_tables"),
  additional_results = file.path(OUTPUT_ROOT, "additional_results"),
  diagnostics = file.path(OUTPUT_ROOT, "diagnostics"),
  metadata = file.path(OUTPUT_ROOT, "metadata")
)

invisible(lapply(c(OUTPUT_ROOT, unname(unlist(OUTPUT_DIRS))), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

RUN_PROFILE <- Sys.getenv("GOLAN_RUN_PROFILE", unset = "publication")
if (!RUN_PROFILE %in% c("publication", "check")) {
  stop("GOLAN_RUN_PROFILE must be 'publication' or 'check'.")
}

SEED <- 42L
N_RANDOMIZATIONS_PART1 <- if (RUN_PROFILE == "publication") 9999L else 3L
N_RANDOMIZATIONS_PART4 <- if (RUN_PROFILE == "publication") 999L else 3L
N_RANDOMIZATIONS_PART6 <- if (RUN_PROFILE == "publication") 9999L else 3L
N_RANDOMIZATIONS_DHAT1 <- if (RUN_PROFILE == "publication") 999L else 3L
N_BOOT_NOCTURNALITY <- if (RUN_PROFILE == "publication") 4999L else 9L
N_BOOT_INTERSPECIES_OVERLAP <- if (RUN_PROFILE == "publication") 499L else 9L
N_BOOT_DHAT1 <- if (RUN_PROFILE == "publication") 1999L else 9L
DELTA1_THRESHOLD <- 50L
N_ACTIVITY_GRID <- 481L
N_Q_EXTREME_CAMERAS <- 15L
requested_cores <- suppressWarnings(as.integer(Sys.getenv("GOLAN_CORES", unset = "")))
available_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(available_cores)) available_cores <- parallel::detectCores()
if (!is.finite(available_cores)) available_cores <- 2L
N_CORES <- if (is.finite(requested_cores) && requested_cores > 0) {
  requested_cores
} else {
  # Four workers are fast for the compiled kernel backend and avoid repeated
  # PSOCK connection exhaustion during low-event permutation tests on Windows.
  max(1L, min(4L, available_cores - 1L))
}

SPECIES <- c(
  "Canis lupus",
  "Canis aureus",
  "Sus scrofa",
  "Gazella gazella",
  "Vulpes vulpes",
  "Meles meles",
  "Hystrix indica"
)

SPECIES_COMMON <- c(
  "Canis lupus" = "Grey wolf",
  "Canis aureus" = "Golden jackal",
  "Sus scrofa" = "Wild boar",
  "Gazella gazella" = "Mountain gazelle",
  "Vulpes vulpes" = "Red fox",
  "Meles meles" = "Caucasian badger",
  "Hystrix indica" = "Indian crested porcupine"
)

SPECIES_COMMON_WRAPPED <- c(
  "Canis lupus" = "Grey wolf",
  "Canis aureus" = "Golden\njackal",
  "Sus scrofa" = "Wild boar",
  "Gazella gazella" = "Mountain\ngazelle",
  "Vulpes vulpes" = "Red fox",
  "Meles meles" = "Caucasian\nbadger",
  "Hystrix indica" = "Indian crested\nporcupine"
)

DETECTION_SPECIES_ORDER <- c(
  "Canis lupus",
  "Canis aureus",
  "Sus scrofa",
  "Gazella gazella",
  "Hystrix indica",
  "Vulpes vulpes",
  "Meles meles"
)

SPECIES_ABBREVIATED <- c(
  "Canis lupus" = "C. lupus",
  "Canis aureus" = "C. aureus",
  "Sus scrofa" = "S. scrofa",
  "Gazella gazella" = "G. gazella",
  "Vulpes vulpes" = "V. vulpes",
  "Meles meles" = "M. meles",
  "Hystrix indica" = "H. indica"
)

LAND_USE_LEVELS <- c(
  "High culling",
  "Low culling",
  "Nature reserves",
  "Minefields"
)

LAND_USE_GRADIENT <- c(
  "High culling",
  "Low culling",
  "Nature reserves",
  "Minefields"
)

LAND_USE_CODES <- c(
  "High culling" = "HC",
  "Low culling" = "LC",
  "Nature reserves" = "NR",
  "Minefields" = "MF"
)

LAND_USE_COLORS <- c(
  "High culling" = "#B85A68",
  "Low culling" = "#3F7FB5",
  "Nature reserves" = "#5E8C31",
  "Minefields" = "#EE7733"
)

DETECTION_FIGURE_COLORS <- c(
  "High culling" = "#B85A68",
  "Low culling" = "#4B8FC0",
  "Nature reserves" = "#8FCB61",
  "Minefields" = "#EE7733"
)

HUMAN_ACTIVITY_COLORS <- c(
  "Low human activity (Q1)" = "#0072B2",
  "High human activity (Q4)" = "#D55E00"
)

RESERVE_COLORS <- c(
  "Temporally zoned reserve" = "#117733",
  "Open-access reserve" = "#D9A214"
)

HUMAN_ALIASES <- c(
  "Bicycle Rider",
  "Vehicle civilian",
  "Homo sapiens-solider",
  "Vehicle construction",
  "Vehicle military",
  "Homo sapiens"
)

PART4_PRESSURES <- c(
  "human activity" = "human.abundance",
  "Wolf activity" = "Wolf.Abundance",
  "jackal activity" = "jackal.abundance",
  "boar activity" = "boar.abundance",
  "Wolf.Culling" = "Wolf.Culling",
  "Jackal.Culling" = "Jackal.Culling",
  "Boar.Culling" = "Boar.Culling"
)

PART4_SPECIES_PRESSURES <- tibble::tribble(
  ~species, ~pressure,
  "Canis lupus", "Wolf.Culling",
  "Canis lupus", "human activity",
  "Canis aureus", "Jackal.Culling",
  "Canis aureus", "human activity",
  "Canis aureus", "Wolf activity",
  "Sus scrofa", "Boar.Culling",
  "Sus scrofa", "human activity",
  "Sus scrofa", "Wolf activity",
  "Gazella gazella", "human activity",
  "Gazella gazella", "Wolf activity",
  "Gazella gazella", "jackal activity",
  "Gazella gazella", "boar activity",
  "Hystrix indica", "human activity",
  "Hystrix indica", "Wolf activity",
  "Hystrix indica", "jackal activity",
  "Hystrix indica", "boar activity",
  "Vulpes vulpes", "human activity",
  "Vulpes vulpes", "Wolf activity",
  "Vulpes vulpes", "jackal activity",
  "Vulpes vulpes", "boar activity",
  "Meles meles", "human activity",
  "Meles meles", "Wolf activity",
  "Meles meles", "jackal activity",
  "Meles meles", "boar activity"
)

# This preserves the final requested badger panel in Figure 3. The validation
# report explicitly notes that it is not the minimum-overlap badger contrast.
FIGURE3_FORCED_PAIR <- c(
  species = "Meles meles",
  group1 = "High culling",
  group2 = "Low culling"
)
