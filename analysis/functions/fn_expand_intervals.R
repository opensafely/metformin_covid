####
## Custom-made function to expand a one-person per row dataset to a event-level dataset
## Taking into account outcome and censoring event dates
## Weekly and monthly intervals (for now). Weekly intervals fixed to 7 days, monthly intervals fixed to 30 days
####
fn_expand_intervals <- function(data, studyend_date, stop_date_columns, outcome_variable, interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)  # Ensure only valid options are chosen
  
  data %>%
    rowwise() %>%
    mutate(
      # Calculate stop_date as the earliest of the defined dates
      stop_date = min(c(!!!syms(stop_date_columns), studyend_date), na.rm = TRUE),
      
      # Define the outcome based on the stop_date and the outcome variable
      outcome = case_when(
        stop_date == .data[[outcome_variable]] ~ 1,  # Outcome 1 for the event defined in outcome_variable
        TRUE ~ 0  # Otherwise, Outcome 0
      )
    ) %>%
    ungroup() %>%
    
    # Create a list of dates based on the interval_type (weekly or monthly)
    mutate(intervals_list = map2(landmark_date, stop_date, ~ {
      if (interval_type == "week") {
        seq(from = .x, to = .y, by = "week")  # Weekly intervals
      } else if (interval_type == "month") {
        seq(from = .x, to = .y, by = "month")  # Monthly intervals
      }
    })) %>%
    unnest(intervals_list) %>%
    
    # Group by patient_id and assign interval-specific information
    group_by(patient_id) %>%
    mutate(
      interval = row_number() - 1,  # Start interval count at 0
      
      # Assign outcome in the final interval (either week or month)
      outcome = if_else(intervals_list <= stop_date & stop_date < intervals_list + ifelse(interval_type == "week", 7, 30), outcome, NA_real_)
    ) %>%
    
    # Fill previous NAs with 0 for consistency in outcomes
    fill(outcome, .direction = "down") %>%
    
    # Replace any remaining NAs with 0 for all rows
    replace_na(list(outcome = 0)) %>%
    
    ungroup() %>%
    
    # Rename the columns dynamically based on the interval type
    rename(
      !!paste0("start_date_", interval_type) := intervals_list,
      !!interval_type := interval  # Assign either "week" or "month" as column names
    )
}