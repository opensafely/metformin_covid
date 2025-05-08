####
# Custom-made function to expand a one-person per row dataset to a event-level dataset

# Follow these rules:
# 1) Expand until first stopping event (Outcome, Censor, Competing) happens and assign intervals (monthly or weekly)
# 2) Assign Outcome/Censor/Competing columns:
## a) If outcome is the stopping event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=1, censor=0, comp_event=0
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

## b) If censoring event is the stopping event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=NA, censor=1, comp_event=NA
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

## c) If competing event is the stopping event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=NA, censor=0, comp_event=1
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

# 3) if outcome and censoring (i.e. LTFU) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 4) if outcome and competing event (i.e. non-covid death) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 5) if censoring and competing event event have the EXACT same date, then pick the competing event as the defining event. Ensured by case_when() order below.

# Important: Since all 3 events above are stopping events, one person-interval can only contain 1 stopping event (outcome, censor, comp)
# Currently, calendar month and calendar week implemented.

####

fn_expand_intervals <- function(data, start_date_variable, stop_date_columns, studyend_date, outcome_date_variable, comp_date_variable, censor_date_variable, interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)  # Ensure only valid options are chosen
  
  data %>%
    rowwise() %>%
    mutate(
      # Assign start date of interval expansion (in our case: landmark_date) & allow for both: import of a date defined separately, or a string name variable part of the input dataset
      start_date = if (is.character(start_date_variable)) .data[[start_date_variable]] else start_date_variable,
      
      # Calculate stop_date as the earliest of the defined dates, including administrative end of study (defined separately)
      stop_date = pmin(!!!syms(stop_date_columns), studyend_date, na.rm = TRUE),

    ) %>%
    ungroup() %>%
    
    # Create a list of dates based on the interval_type (weekly or monthly)
    mutate(intervals_list = map2(start_date, stop_date, ~ {
      last_possible_date <- .y - days(1) # Subtract one day; left-closed, right-open
      if (interval_type == "week") {
        seq(from = .x, to = last_possible_date, by = "week")
      } else {
        seq(from = .x, to = last_possible_date, by = "month")
      }
    })) %>%
    unnest(intervals_list) %>%
    
    # Group by patient_id and assign interval-specific information
    group_by(patient_id) %>%
    
    # End date interval
    mutate(
      !!paste0("end_date_", interval_type) := lead(intervals_list, default = stop_date[1] + 1) - days(1)
    ) %>%
    
    mutate(
      interval = row_number() - 1, # Start interval count at 0
      
      # Flag end event interval (i.e. where outcome, competing, or censoring occurs)
      is_event_interval = intervals_list <= stop_date & stop_date <= (!!sym(paste0("end_date_", interval_type))),  
      
      # Flag what kind of end event ocurred
      is_outcome_event = if_else(is_event_interval & stop_date == .data[[outcome_date_variable]], 1L, 0L),
      is_cens_event = if_else(is_event_interval & stop_date == .data[[censor_date_variable]], 1L, 0L),
      is_comp_event = if_else(is_event_interval & stop_date == .data[[comp_date_variable]], 1L, 0L),
      # note: if the stopping event is the end of follow-up, then NA is assigned to the corresponding person-interval
      
      # Assign outcome column, see rules above
      outcome = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 1, # Is the outcome stopping event happening in next interval? If so, assign 1 to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ NA_real_, # Is the censoring stopping event happening in next interval? If so, assign NA to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ NA_real_, # Is the competing stopping event happening in next interval? If so, assign NA to current (= prior) interval
        # define entry for interval when event is happening
        is_outcome_event == 1 ~ NA_real_, # to current interval assign NA (i.e. to interval where death happens)
        is_event_interval & stop_date == .data[[comp_date_variable]] ~ NA_real_, # stopping event is the competing event, assign NA to current interval
        is_event_interval & stop_date == .data[[censor_date_variable]] ~ NA_real_, # stopping event is the censoring event, assign NA to current interval
        is_event_interval & stop_date == studyend_date ~ 0, # stopping event is end of follow-up, assign 0 to current interval
        TRUE ~ 0 # Assign 0 to all other intervals
        ),
      # Assign competing event column, see rules above
      comp_event = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 0, # Is the outcome stopping event happening in next interval? If so, assign 0 to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ NA_real_, # Is the censoring stopping event happening in next interval? If so, assign NA to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ 1, # Is the competing stopping event happening in next interval? If so, assign 1 to current (= prior) interval
        # define entry for interval when event is happening
        is_comp_event == 1 ~ NA_real_,
        is_event_interval & stop_date == .data[[outcome_date_variable]] ~ NA_real_,
        is_event_interval & stop_date == .data[[censor_date_variable]] ~ NA_real_,
        is_event_interval & stop_date == studyend_date ~ 0,
        TRUE ~ 0
        ),
      # Assign censor event column, see rules above
      censor = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 0, # Is the outcome stopping event happening in next interval? If so, assign 0 to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ 1, # Is the censoring stopping event happening in next interval? If so, assign 1 to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ 0, # Is the competing stopping event happening in next interval? If so, assign 0 to current (= prior) interval
        # define entry for interval when event is happening
        is_cens_event == 1 ~ NA_real_,
        is_event_interval & stop_date == .data[[outcome_date_variable]] ~ NA_real_,
        is_event_interval & stop_date == .data[[comp_date_variable]] ~ NA_real_,
        is_event_interval & stop_date == studyend_date ~ 0,
        TRUE ~ 0
      ),
      # check if anyone has more than 1 stopping event interval
      followup_stop = cumsum(is_event_interval) > 1
    ) %>%

    ungroup() %>%
    
    # Drop rows where outcome, censor and comp_event are all NA -> these are the intervals where outcomes happened, but we shifted the info to the row above => drop these rows
    filter(!(is.na(outcome) & is.na(censor) & is.na(comp_event))) %>%
    
    # Rename the columns dynamically based on the interval type
    rename(
      !!paste0("start_date_", interval_type) := intervals_list,
      !!interval_type := interval  # Assign either "week" or "month" as column names
    )
}