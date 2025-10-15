####
# Custom-made function to expand a one-person per row dataset to a event-level dataset
# Expand until first stopping event happens and assign calendar intervals (monthly or weekly)
# Required packages (load in main script):
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(rlang)
# library(lubridate)
# Or: library(tidyverse) & library(lubridate) & library(rlang) 

fn_expand_intervals <- function(data, start_date_variable, stop_date_columns, studyend_date, interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)  # Ensure only valid options are chosen
  
  data %>%
    rowwise() %>%
    mutate(
      start_date = if (is.character(start_date_variable)) .data[[start_date_variable]] else start_date_variable,
      stop_date = pmin(!!!syms(stop_date_columns), studyend_date, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    
    # Create list of (start, stop) pairs for each interval
    mutate(intervals_list = map2(start_date, stop_date, ~ {
      current_start <- .x
      intervals <- list()
      while (current_start <= .y) {
        next_start <- if (interval_type == "week") current_start + weeks(1) else current_start %m+% months(1)
        current_stop <- next_start - days(1)
        if (current_stop > .y) current_stop <- .y
        intervals <- append(intervals, list(list(start = current_start, stop = current_stop)))
        current_start <- next_start
      }
      intervals
    })) %>%
    
    # Unnest the list into rows and widen start/stop
    unnest(intervals_list) %>%
    unnest_wider(intervals_list, names_sep = "_") %>%
    
    group_by(patient_id) %>%
    mutate(
      interval = row_number() - 1 # Start interval count at 0
    ) %>%
    ungroup() %>%
    
    # Rename the interval columns
    rename(
      !!paste0("start_date_", interval_type) := intervals_list_start,
      !!paste0("end_date_", interval_type) := intervals_list_stop,
      !!interval_type := interval
    )
}
