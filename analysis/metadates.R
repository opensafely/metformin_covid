# # # # # # # # # # # # # # # # # # # # #
# This script creates the study dates (json file) to be called via a yaml action
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library('tidyverse')
library('here')

# create study_dates ----
study_dates <-
  list(
    studystart_date = "2020-01-01", # start of pandemic
    studyend_date = "2022-04-01" # end of mass testing
  )

jsonlite::write_json(study_dates, path = "output/study_dates.json", auto_unbox = TRUE, pretty=TRUE)