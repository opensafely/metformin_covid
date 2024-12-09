# # # # # # # # # # # # # # # # # # # # #
# This script creates the study dates (json file) to be called via a yaml action
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library('tidyverse')
library('here')

# create study_dates ----
study_dates <-
  list(
    studystart_date = "2020-02-01", #start of pandemic
    studyend_date = "2022-04-01", #end of mass testing
    feasibilityend_date = "2023-12-31", #end date of feasibility study
    studystartlandmark_date = "2018-07-01", #start date of landmark study
    landmark_date = "2020-02-01", #landmark = pandemic start = follow-up start
    mid2018_date = "2018-07-01"
  )

jsonlite::write_json(study_dates, path = "output/study_dates.json", auto_unbox = TRUE, pretty=TRUE)