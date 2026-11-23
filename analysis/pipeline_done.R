################################################################################
# Collector action: writes a marker file once all actions of a pipeline ran.
# Run a whole pipeline via its collector action, e.g.:
#   opensafely run run_landmark (original landmark pipeline)
#   opensafely run st_run_all (sequential trials pipeline)
# Usage: r:v2 analysis/pipeline_done.R <landmark|seqtrials>
################################################################################
args <- commandArgs(trailingOnly = TRUE)
pipeline <- if (length(args) == 0) "landmark" else args[[1]]

out_file <- switch(pipeline,
  landmark = here::here("output", "landmark_done.txt"),
  seqtrials = here::here("output", "seqtrials", "st_done.txt"),
  stop("Unknown pipeline: '", pipeline, "'. Use 'landmark' or 'seqtrials'.")
)

fs::dir_create(dirname(out_file))
writeLines(paste0("Pipeline '", pipeline, "' completed."), out_file)
