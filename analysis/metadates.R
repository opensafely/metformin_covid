################################################################################
# This script creates the study dates (json file) to be called via a yaml action
# Usage:
#   r:v2 analysis/metadates.R -> landmark pipeline:  output/study_dates.json
#   r:v2 analysis/metadates.R seqtrials -> seqtrials pipeline: output/seqtrials/study_dates.json
################################################################################

################################################################################
# 0.0 Import libraries
################################################################################
library('jsonlite')
library('here')

################################################################################
# 0.1 Pipeline argument (default: landmark, i.e. the original pipeline)
################################################################################
# This script is also source()d by many analysis scripts to get `study_dates`.
# Only read the command line when run directly as an action (Rscript), otherwise
# commandArgs() would return the calling script's arguments (e.g. "no_interaction").
# When sourced, the pipeline defaults to landmark; a calling script can set
# `metadates_pipeline <- "seqtrials"` before source() to get the seqtrials dates.
# Everything is wrapped in local() so that only `study_dates` is created in the caller.
.metadates_sourced <- sys.nframe() > 0L # 0 when run directly via Rscript, > 0 via source()
study_dates <- local({
  if (!.metadates_sourced) {
    args <- commandArgs(trailingOnly = TRUE)
    pipeline <- if (length(args) == 0) "landmark" else args[[1]]
  } else {
    pipeline <- get0("metadates_pipeline", ifnotfound = "landmark")
  }

  ##############################################################################
  # 1 Define study_dates per pipeline and write as json to output
  ##############################################################################
  if (pipeline == "landmark") {
    out_dir <- here::here("output")
    study_dates <-
      list(
        studyend_date = "2022-04-01",
        pandemicstart_date = "2020-02-01",
        mid2018_date = "2018-08-01",
        mid2019_date = "2019-08-01"
      )
  } else if (pipeline == "seqtrials") {
    out_dir <- here::here("output", "seqtrials")
    study_dates <-
      list(
        studystart_date = "2020-02-01", # start of pandemic = first trial (T0)
        studyend_date = "2022-03-01" # end of systematic SARS-CoV-2 mass testing in the UK
      )
  } else {
    stop("Unknown pipeline: '", pipeline, "'. Use 'landmark' (default) or 'seqtrials'.")
  }

  fs::dir_create(out_dir)
  write_json(study_dates, path = file.path(out_dir, "study_dates.json"), auto_unbox = TRUE, pretty = TRUE)
  study_dates
})
rm(.metadates_sourced)
