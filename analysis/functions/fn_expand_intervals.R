####
# Custom-made function to expand a one-person per row dataset to a event-level dataset
# Expand the dataset into intervals and follow these rules:
# a) If outcome (out_date_severecovid_afterlandmark) is reached first, assign outcome=1, censor=0, comp_event=0 to the interval when it happened and stop expanding
# b) If competing event (out_date_death_afterlandmark) is reached first, assign outcome=NA, censor=0, comp_event=1 to the interval when it happened and stop expanding
# c) If censoring event (out_date_ltfu_afterlandmark) is reached first, assign outcome=NA, censor=1, comp_event=NA to the interval when it happened and stop expanding
# d) If studyend_date is reached first, then assign outcome=0, censor=0, comp_event=0 to the interval when it happened and stop expanding
# Weekly and monthly intervals (for now). Weekly intervals fixed to 7 days. Monthly intervals fixed to 31 days (to ensure max possible length of a month)
####
fn_expand_intervals <- function(data, start_date_variable, stop_date_columns, studyend_date, outcome_date_variable, comp_date_variable, censor_date_variable, interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)  # Ensure only valid options are chosen
  
  data %>%
    rowwise() %>%
    mutate(
      # Assign start date of interval expansion (in our case: landmark_date)
      start_date = .data[[start_date_variable]],
      
      # Calculate stop_date as the earliest of the defined dates, including administrative end of study (usually defined separately)
      stop_date = pmin(!!!syms(stop_date_columns), studyend_date, na.rm = TRUE),
      
      # Define the outcome variable based on which outcome date variable defined the stop_date 
      outcome = case_when(
        stop_date == .data[[outcome_date_variable]] ~ 1,
        stop_date == .data[[comp_date_variable]] ~ NA_real_,
        stop_date == .data[[censor_date_variable]] ~ NA_real_,
        stop_date == studyend_date ~ 0,
        TRUE ~ NA_real_ # but is illogical, can't happen
      ),
      
      # Define the competing event variable based on which outcome date variable defined the stop_date
      comp_event = case_when(
        stop_date == .data[[outcome_date_variable]] ~ 0,
        stop_date == .data[[comp_date_variable]] ~ 1,
        stop_date == .data[[censor_date_variable]] ~ NA_real_,
        stop_date == studyend_date ~ 0,
        TRUE ~ NA_real_ # but is illogical, can't happen
      ),
      
      # Define the censoring event variable based on which outcome date variable defined the stop_date
      censor = case_when(
        stop_date == .data[[outcome_date_variable]] ~ 0,
        stop_date == .data[[comp_date_variable]] ~ 0,
        stop_date == .data[[censor_date_variable]] ~ 1,
        stop_date == studyend_date ~ 0,
        TRUE ~ NA_real_ # but is illogical, can't happen
      )
      
    ) %>%
    ungroup() %>%
    
    # Create a list of dates based on the interval_type (weekly or monthly)
    mutate(intervals_list = map2(start_date, stop_date, ~ {
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
      # Assign outcome variable in the final interval
      outcome = if_else(intervals_list <= stop_date & stop_date < intervals_list + ifelse(interval_type == "week", 7, 31), outcome, NA_real_),
      # Assign competing event variable in the final interval
      comp_event = if_else(intervals_list <= stop_date & stop_date < intervals_list + ifelse(interval_type == "week", 7, 31), comp_event, NA_real_),
      # Assign cesnoring event variable in the final interval
      censor = if_else(intervals_list <= stop_date & stop_date < intervals_list + ifelse(interval_type == "week", 7, 31), censor, NA_real_)
    ) %>%
    
    # Fill previous NAs with 0
    mutate(
      outcome = if_else(is.na(outcome) & interval < max(interval), 0, outcome),
      comp_event = if_else(is.na(comp_event) & interval < max(interval), 0, comp_event),
      censor = if_else(is.na(censor) & interval < max(interval), 0, censor)
    ) %>%

    ungroup() %>%
    
    # Rename the columns dynamically based on the interval type
    rename(
      !!paste0("start_date_", interval_type) := intervals_list,
      !!interval_type := interval  # Assign either "week" or "month" as column names
    )
}
