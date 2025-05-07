####
# Custom-made function to expand a one-person per row dataset to a event-level dataset

# Expand the dataset into intervals and follow these rules:
# a) If outcome (out_date_severecovid_afterlandmark) is reached first, stop expanding, define the interval the event happened (but ignore that interval), and assign to the PREVIOUS interval: outcome=1, censor=0, comp_event=0
# => this becomes an indicator for the outcome event in interval𝑘 + 1! See Boston CAUSALab material, they used the same rules. See 'vac_random' dataset in first course session.
# b) If competing event (out_date_death_afterlandmark) is reached first, stop expanding, define the interval the event happened and assign outcome=NA, censor=0, comp_event=1 to that interval
# c) If censoring event (out_date_ltfu_afterlandmark) is reached first, stop expanding, define the interval the event happened and assign outcome=NA, censor=1, comp_event=NA to that interval
# d) If studyend_date is reached first, then assign outcome=0, censor=0, comp_event=0 to the interval when it happened and stop expanding

# further rules:
# e) if outcome and censoring (i.e. LTFU) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# f) if outcome and competing event (i.e. non-covid death) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# g) if censoring and competing event event have the EXACT same date, then pick the competing event as the defining event. Ensured by case_when() order below.

# Weekly and monthly intervals (for now). Weekly intervals fixed to 7 days. Monthly intervals fixed to 31 days (to ensure max possible length of a month)
####
fn_expand_intervals <- function(data, start_date_variable, stop_date_columns, studyend_date, outcome_date_variable, comp_date_variable, censor_date_variable, interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)  # Ensure only valid options are chosen
  
  data %>%
    rowwise() %>%
    mutate(
      # Assign start date of interval expansion (in our case: landmark_date & allow for both: direct import of a date, or a string name of the dataset)
      start_date = if (is.character(start_date_variable)) .data[[start_date_variable]] else start_date_variable,
      
      # Calculate stop_date as the earliest of the defined dates, including administrative end of study (usually defined separately)
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
    
    # Rename/explicitly show intervals
    mutate(
      !!paste0("end_date_", interval_type) := lead(intervals_list, default = stop_date[1] + 1) - days(1)
    ) %>%
    
    mutate(
      interval = row_number() - 1, # Start interval count at 0
      
      # Flag end event interval (i.e. where outcome, competing, or censoring occurs)
      is_event_interval = intervals_list <= stop_date & stop_date <= (!!sym(paste0("end_date_", interval_type))),  # left-closed, right-open
      
      # Flag what kind of end event ocurred
      is_outcome_event = if_else(is_event_interval & stop_date == .data[[outcome_date_variable]], 1L, 0L),
      is_cens_event = if_else(is_event_interval & stop_date == .data[[censor_date_variable]], 1L, 0L),
      is_comp_event = if_else(is_event_interval & stop_date == .data[[comp_date_variable]], 1L, 0L),
      
      # Assign outcome column, see rules a-d above
      outcome = case_when(
        lead(is_outcome_event, default = 0) == 1 ~ 1, # Is event happening in next interval? If so, assign 1
        is_outcome_event == 1 ~ NA_real_, # to current interval assign NA
        is_event_interval & stop_date == .data[[comp_date_variable]] ~ NA_real_, # stopping event is the competing event, assign NA to 'outcome' to current interval, see rules above
        is_event_interval & stop_date == .data[[censor_date_variable]] ~ NA_real_, # stopping event is the censoring event, assign NA to 'outcome' to current interval, see rules above
        is_event_interval & stop_date == studyend_date ~ 0, # stopping event is end of follow-up, assign 0 to 'outcome' to current interval, see rules above
        TRUE ~ 0
        ),
      # Assign competing event column, see rules a-d above
      comp_event = case_when(is_event_interval & stop_date == .data[[outcome_date_variable]] ~ NA_real_, # stopping event is an outcome, assign NA to comp_event to current interval (as per rule a) and previous interval will automatically have O, see default
                          is_event_interval & stop_date == .data[[comp_date_variable]] ~ 1, # stopping event is the competing event, assign 1 to 'comp_event' to current interval
                          is_event_interval & stop_date == .data[[censor_date_variable]] ~ NA_real_, # stopping event is an censoring event, assign NA to 'comp_event' to current interval
                          is_event_interval & stop_date == studyend_date ~ 0, # stopping event is end of follow-up
                          TRUE ~ 0 # Assign 0 to all other intervals
                          ),
      # Assign censor event column, see rules a-d above
      censor = case_when(is_event_interval & stop_date == .data[[outcome_date_variable]] ~ NA_real_, # stopping event is an outcome, assign NA to comp_event to current interval (as per rule a) and previous interval will automatically have O, see defaul
                             is_event_interval & stop_date == .data[[comp_date_variable]] ~ 0, # stopping event is a competing event, assign 0 to 'censor' to current interval
                             is_event_interval & stop_date == .data[[censor_date_variable]] ~ 1, # stopping event is the censoring event, assign 1 to 'censor' to current interval
                             is_event_interval & stop_date == studyend_date ~ 0, # stopping event is end of follow-up
                             TRUE ~ 0 # Assign 0 to all other intervals
                             ),
      followup_stop = cumsum(is_event_interval) > 1
    ) %>%
    
    # If e.g. stop_date is 2022-04-01 and start_date_month is 2022-03-02 it will start a new unnecessary interval, hence is_event_interval will be assigned to the interval it happened and to a next unnecessary interval. Drop that next one.
    filter(!followup_stop) %>%
    ungroup() %>%
    
    # Drop rows where outcome, censor and comp_event are all NA -> these are the intervals where outcomes happened, but we shifted the info to the row above => drop this row.
    filter(!(is.na(outcome) & is.na(censor) & is.na(comp_event))) %>%
    # this is an edge case, whereby is_outcome_interval was assigned to two rows due to stop_date is 2022-04-01 and start_date_month is 2022-03-02. 
    # This will lead to: a) a row with outcome==1, censor==0, comp_event==0 (the row we want to keep, using info from k + 1)
    # b) a row with outcome==NA, censor==NA, comp_event==NA -> excluded above
    # c) a row with outcome==1, censor==NA, comp_event==NA -> exclude now
    filter(!(is.na(censor) & is.na(comp_event))) %>% 
    
    # Rename the columns dynamically based on the interval type
    rename(
      !!paste0("start_date_", interval_type) := intervals_list,
      !!interval_type := interval  # Assign either "week" or "month" as column names
    )
}
