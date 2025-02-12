################################################################################
## This script does only the following:
# If the OpenSAFELY is run locally, it directs it to data_process_dummydata.R
# otherwise (if run on the real data), it directs it to data_process.R
################################################################################

if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  source("analysis/data_process_dummydata.R")
  message("Dummy data pathway")
} else {
  source("analysis/data_process.R")
  message("Real data pathway")
}