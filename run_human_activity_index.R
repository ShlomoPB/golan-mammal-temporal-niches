source(file.path("R", "00_config.R"))
source(file.path("R", "human_activity_index.R"))

args <- commandArgs(trailingOnly = TRUE)
valid_args <- c("--refit-model", "--refit-candidates")
unknown_args <- setdiff(args, valid_args)
if (length(unknown_args) > 0L) {
  stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "))
}

refit_candidates <- "--refit-candidates" %in% args
refit_model <- refit_candidates || "--refit-model" %in% args
message(
  "Human-activity reconstruction; refit model: ",
  refit_model,
  "; refit candidates: ",
  refit_candidates
)
run_human_activity_index(
  refit_model = refit_model,
  refit_candidates = refit_candidates
)
