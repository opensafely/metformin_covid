# Load packages ----------------------------------------------------------------
print('Load packages')

library(magrittr)

# Source common functions ------------------------------------------------------
print('Source common functions')

source(here::here("analysis", "functions", "utility.R"))

# Define model output folder ---------------------------------------
print("Creating output/model output folder")

# setting up the sub directory
makeout_dir <- "output/make_output/"
model_dir <- "output/"

# check if sub directory exists, create if not
fs::dir_create(here::here(makeout_dir))

# List available model outputs -------------------------------------------------
print('List available model outputs')

files_all <- list.files(model_dir, pattern = "^results_cox_.*\\.csv$", full.names = FALSE)
files_R <- files_all[!grepl("_ra|_minset", files_all)] # Filter out any filenames that contain "_ra" or "_minset" from earlier naming structure

# Combine R model output -------------------------------------------------------
print('Combine R model output')

df <- NULL

for (i in files_R) {
  ## Load model output
  tmp <- readr::read_csv(paste0(model_dir, i))

  ## Handle errors
  if (colnames(tmp)[1] == "error") {
    dummy <- data.frame(
      model = "",
      exposure = "",
      outcome = gsub(".*-", "", gsub(".csv", "", i)),
      term = "",
      lnhr = NA,
      se_lnhr = NA,
      hr = NA,
      conf_low = NA,
      conf_high = NA,
      N_total = NA,
      N_exposed = NA,
      N_events = NA,
      person_time_total = NA,
      outcome_time_median = NA,
      strata_warning = "",
      surv_formula = "",
      input = "",
      error = tmp$error
    )
    tmp <- dummy
  } else {
    tmp$error <- ""
  }

  ## Add source file name
  tmp$name <- gsub("results_cox_", "", gsub(".csv", "", i))

  ## Append to master dataframe
  df <- plyr::rbind.fill(df, tmp)
}

# Save model output ------------------------------------------------------------
print('Save model output')

df <- df[, c(
  "name",
  "error",
  "model",
  "term",
  "lnhr",
  "se_lnhr",
  "hr",
  "conf_low",
  "conf_high",
  "N_total",
  "N_exposed",
  "N_events",
  "person_time_total",
  "outcome_time_median",
  "strata_warning",
  "surv_formula"
)]

readr::write_csv(df, paste0(makeout_dir, "results_cox.csv"))

# Perform redaction ------------------------------------------------------------
print('Perform redaction')

df$N_total_midpoint6 <- fn_roundmid_any(df$N_total)
df$N_exposed_midpoint6 <- fn_roundmid_any(df$N_exposed)
df$N_events_midpoint6 <- fn_roundmid_any(df$N_events)
df[, c("N_total", "N_exposed", "N_events")] <- NULL

# Save model output ------------------------------------------------------------
print('Save model output')

readr::write_csv(
  df,
  paste0(makeout_dir, "results_cox_midpoint6.csv")
)
