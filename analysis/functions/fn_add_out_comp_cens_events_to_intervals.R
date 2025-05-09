####
# Custom-made function to add outcome, competing and censoring events to an expanded person-interval dataset based on fn_expand_intervals

# Follow these rules:
# Assign Outcome/Censor/Competing columns:
## a) If outcome is the event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=1, censor=0, comp_event=0
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

## b) If censoring event is the event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=NA, censor=1, comp_event=NA
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

## c) If competing event is the event: 
### i) assign to the PREVIOUS interval row/person-interval: outcome=NA, censor=0, comp_event=1
### ii) assign to CURRENT interval row/person-interval: outcome=NA, censor=NA, comp_event=NA

# 3) if outcome and censoring (i.e. LTFU) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 4) if outcome and competing event (i.e. non-covid death) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 5) if censoring and competing event event have the EXACT same date, then pick the competing event as the defining event. Ensured by case_when() order below.

# Important: 1 person-interval can only contain 1 event (outcome, censor, comp)
# Currently, calendar month and calendar week implemented.

####
fn_add_out_comp_cens_events_to_intervals <- function(data, outcome_date_variable, comp_date_variable, censor_date_variable, studyend_date, interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)
  
  data %>%
    group_by(patient_id) %>%
    mutate(
      event_date = pmin(.data[[outcome_date_variable]], .data[[comp_date_variable]], .data[[censor_date_variable]], studyend_date, na.rm = TRUE),
      
      is_event_interval = (.data[[paste0("start_date_", interval_type)]] <= event_date) & (event_date <= .data[[paste0("end_date_", interval_type)]]),
      
      is_outcome_event = if_else(is_event_interval & event_date == .data[[outcome_date_variable]], 1L, 0L),
      is_cens_event = if_else(is_event_interval & event_date == .data[[censor_date_variable]], 1L, 0L),
      is_comp_event = if_else(is_event_interval & event_date == .data[[comp_date_variable]], 1L, 0L),
      
      outcome = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 1, # Is the outcome event happening in next interval? If so, assign 1 to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ NA_real_, # Is the censoring event happening in next interval? If so, assign NA to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ NA_real_, # Is the competing event happening in next interval? If so, assign NA to current (= prior) interval
        # define entry for interval when event is happening
        is_outcome_event == 1 ~ NA_real_, # to current interval assign NA (i.e. to interval where death happens)
        is_event_interval & event_date == .data[[comp_date_variable]] ~ NA_real_, # event is the competing event, assign NA to current interval
        is_event_interval & event_date == .data[[censor_date_variable]] ~ NA_real_, # event is the censoring event, assign NA to current interval
        is_event_interval & event_date == studyend_date ~ 0, # event is end of follow-up, assign 0 to current interval
        TRUE ~ 0 # Assign 0 to all other intervals
      ),
      
      comp_event = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 0, # Is the outcome event happening in next interval? If so, assign 0 to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ NA_real_, # Is the censoring event happening in next interval? If so, assign NA to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ 1, # Is the competing event happening in next interval? If so, assign 1 to current (= prior) interval
        # define entry for interval when event is happening
        is_comp_event == 1 ~ NA_real_,
        is_event_interval & event_date == .data[[outcome_date_variable]] ~ NA_real_,
        is_event_interval & event_date == .data[[censor_date_variable]] ~ NA_real_,
        is_event_interval & event_date == studyend_date ~ 0,
        TRUE ~ 0
      ),
      
      censor = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 0, # Is the outcome event happening in next interval? If so, assign 0 to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ 1, # Is the censoring event happening in next interval? If so, assign 1 to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ 0, # Is the competing event happening in next interval? If so, assign 0 to current (= prior) interval
        # define entry for interval when event is happening
        is_cens_event == 1 ~ NA_real_,
        is_event_interval & event_date == .data[[outcome_date_variable]] ~ NA_real_,
        is_event_interval & event_date == .data[[comp_date_variable]] ~ NA_real_,
        is_event_interval & event_date == studyend_date ~ 0,
        TRUE ~ 0
      ),
      # check if anyone has more than 1 event interval
      dupl_interval = cumsum(is_event_interval) > 1
    ) %>%
    
    ungroup() %>%
    # Drop rows where outcome, censor and comp_event are all NA -> these are the intervals where outcomes happened, but we shifted the info to the row above => drop these rows
    filter(!(is.na(outcome) & is.na(censor) & is.na(comp_event)))
}