################################################################################
# A custom made function to add columns period_week and period_month
# Marking the study population in weekly/monthly intervals (1/0) based on covid_test_positive_date (baseline_date)
# Based on and Credits to https://github.com/opensafely/pax-non-users/tree/2dbf044472efdcfeb86f8fc2c8eea222e7eefe32
################################################################################

# Arguments:
# data_extracted: data_frame with study population
# study_dates: study_dates from lib/design/study-dates.json

add_period_cuts <- function(data_extracted, study_dates){
  seq_dates_start_interval_week <-
    seq(study_dates$studystart_date, study_dates$studyend_date + 1, by = "1 week") # adding +1 to the end date ensures that the end date is included in the sequence
  seq_dates_start_interval_month <-
    seq(study_dates$studystart_date, study_dates$studyend_date + 1, by = "1 month")
  seq_dates_start_interval_2month <-
    seq(study_dates$studystart_date, study_dates$studyend_date + 1, by = "2 months")
  seq_dates_start_interval_3month <-
    seq(study_dates$studystart_date, study_dates$studyend_date + 1, by = "3 months")

  data_extracted <- data_extracted %>%
    mutate(period_week = cut(baseline_date,
                             breaks = seq_dates_start_interval_week,
                             include.lowest = TRUE, # if a baseline_date falls exactly on the start date of the interval, it will be included in that interval
                             right = FALSE, # intervals are closed on the left side and open on the right side
                             labels = 1:117), # adapt according to study period! Currently 27 months / 116 weeks (1.1.2020-1.4.2022)
           period_month = cut(baseline_date,
                              breaks = seq_dates_start_interval_month,
                              include.lowest = TRUE,
                              right = FALSE,
                              labels = 1:27),
           period_2month = cut(baseline_date,
                              breaks = seq_dates_start_interval_2month,
                              include.lowest = TRUE,
                              right = FALSE,
                              labels = 1:13.5),
           period_3month = cut(baseline_date,
                               breaks = seq_dates_start_interval_3month,
                               include.lowest = TRUE,
                               right = FALSE,
                               labels = 1:9))
}
