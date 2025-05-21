####
# Custom-made function to add outcome, competing and censoring events to an expanded person-interval dataset
# CAVE: use it only after the dataset has been interval-expanded (fn_expand_intervals.R) and after all one-time events have been added (fn_add_firstever_events_to_intervals.R)

# Assign columns called "outcome", "censor", "comp_event" as following:

## a) If outcome is the row/person-interval event: 
### i) assign to the PREVIOUS row/person-interval: outcome=1, censor=0, comp_event=0
### ii) assign to CURRENT row/person-interval: outcome=NA, censor=NA, comp_event=NA

## b) If censoring event is the row/person-interval event: 
### i) assign to the PREVIOUS row/person-interval: outcome=NA, censor=1, comp_event=NA
### ii) assign to CURRENT row/person-interval: outcome=NA, censor=NA, comp_event=NA

## c) If competing event is the row/person-interval event: 
### i) assign to the PREVIOUS row/person-interval: outcome=NA, censor=0, comp_event=1
### ii) assign to CURRENT row/person-interval: outcome=NA, censor=NA, comp_event=NA

# 3) if outcome and censoring (i.e. LTFU) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 4) if outcome and competing event (i.e. non-covid death) event have the EXACT same date, then pick the outcome as the defining event. Ensured by case_when() order below.
# 5) if censoring and competing event event have the EXACT same date, then pick the competing event as the defining event. Ensured by case_when() order below.

# 6) Delete all rows/person-intervals that happen after the outcome

# Currently, only calendar month and calendar week are implemented

fn_add_and_shift_out_comp_cens_events <- function(data, 
                                          outcome_date_variable, 
                                          comp_date_variable, 
                                          censor_date_variable, 
                                          studyend_date, 
                                          interval_type = c("week", "month")) {
  interval_type <- match.arg(interval_type)
  
  data %>%
    group_by(patient_id) %>%
    mutate(
      # I need to use the dates (not the bin flags), since in some scenarios an outcome and censoring or competing event might occur in same interval. I will only consider what happens first (pmin)
      event_date = pmin(.data[[outcome_date_variable]], .data[[comp_date_variable]], .data[[censor_date_variable]], studyend_date, na.rm = TRUE),

      is_event_interval = (.data[[paste0("start_date_", interval_type)]] <= event_date) & (event_date <= .data[[paste0("end_date_", interval_type)]]),

      is_outcome_event = if_else(is_event_interval & event_date == .data[[outcome_date_variable]], 1L, 0L),
      is_comp_event = if_else(is_event_interval & event_date == .data[[comp_date_variable]], 1L, 0L),
      is_cens_event = if_else(is_event_interval & event_date == .data[[censor_date_variable]], 1L, 0L),
      
      outcome = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 1, # Is the outcome event happening in next interval? If so, assign 1 to current (= prior) interval
        lead(is_comp_event, default = 0) == 1 ~ NA_real_, # Is the competing event happening in next interval? If so, assign NA to current (= prior) interval
        lead(is_cens_event, default = 0) == 1 ~ NA_real_, # Is the censoring event happening in next interval? If so, assign NA to current (= prior) interval
        # define entry for interval when event is happening
        is_outcome_event == 1 ~ NA_real_, # to current interval assign NA (i.e. to interval where death happens)
        is_comp_event == 1 ~ NA_real_, # event is the competing event, assign NA to current interval
        is_cens_event == 1 ~ NA_real_, # event is the censoring event, assign NA to current interval
        TRUE ~ 0 # Assign 0 to all other intervals - and in case of no event happened (i.e. end of follow-up)
      ),
      
      comp_event = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 0, 
        lead(is_comp_event, default = 0) == 1 ~ NA_real_, 
        lead(is_cens_event, default = 0) == 1 ~ 1,
        # define entry for interval when event is happening
        is_comp_event == 1 ~ NA_real_,
        is_outcome_event == 1 ~ NA_real_,
        is_cens_event == 1 ~ NA_real_,
        TRUE ~ 0
      ),
      
      censor = case_when(
        # define entry for interval prior to event
        lead(is_outcome_event, default = 0) == 1 ~ 0, 
        lead(is_comp_event, default = 0) == 1 ~ 0,
        lead(is_cens_event, default = 0) == 1 ~ 1, 
        # define entry for interval when event is happening
        is_comp_event == 1 ~ NA_real_,
        is_outcome_event == 1 ~ NA_real_,
        is_cens_event == 1 ~ NA_real_,
        TRUE ~ 0
      )
      
    ) %>%
    
    # In scenarios where the outcome is not a final/fatal outcome (e.g. covid hospitalization) the dataset will have rows after such an event.
    # Make sure there is no data after outcome, comp, or censoring event happened
    mutate(row_num = row_number()) %>%
    mutate(first_event_row = min(ifelse(is_event_interval, row_num, Inf))) %>% # in case there is no event at all (i.e. end of follow-up is the stopping event)
    filter(row_num <= first_event_row) %>%
    
    ungroup() %>%
    
    # Drop rows where outcome, censor and comp_event are all NA -> these are the intervals where outcomes happened, but we shifted the info to the row above => drop these rows
    filter(!(is.na(outcome) & is.na(censor) & is.na(comp_event))) %>%
    
    # Drop helper variables
    select(-event_date, -is_event_interval, -is_outcome_event, -is_comp_event, -is_cens_event, -row_num, -first_event_row)
  
}