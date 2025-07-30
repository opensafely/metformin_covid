####
# Custom-made function to add time-fixed covariates, i.e., those occurring once and remaining stable (these are usually the highest number of variables in a dataset), to person-time interval data
## RULES:
# a) if event date is not NA and happened before the minimum start date of all intervals of a person, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happened after the maximum end date of all intervals of a person, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals (assuming no documentation = no event)
# d) if event date happened during follow-up, then assign 1 (in flag variable) to next person-interval (lag), and 0 to all person-intervals before, and 1 to all person-intervals after (stable/time-fixed event)

fn_assign_time_fixed_cov <- function(df, date_vars, start_var, end_var) {
  
  df_out <- df
  
  walk(date_vars, function(var) {
    # create binary flagging variable by replacing the _date_ with a corresponding _bin_
    bin_var <- str_replace(var, "_date_", "_bin_")
    
    df_out <<- df_out %>%
      group_by(patient_id) %>%
      
      mutate(
        elig_date = {
          # to be on the safe side if the entire variable is NA
          na_value <- if (is.numeric(.data[[var]])) {
            NA_real_
          } else if (inherits(.data[[var]], "Date")) {
            as.Date(NA)
          } else {
            NA
          }
          if (all(is.na(.data[[var]]))) na_value else first(na.omit(.data[[var]]))
        },
        !!bin_var := case_when(
          is.na(elig_date) ~ 0,
          elig_date < min(.data[[start_var]], na.rm = TRUE) ~ 1, # before follow-up
          elig_date >= min(.data[[start_var]], na.rm = TRUE) & elig_date <= max(.data[[end_var]], na.rm = TRUE) ~ NA_real_, # during follow-up; handled in next mutate
          elig_date > max(.data[[end_var]], na.rm = TRUE) ~ 0 # after follow-up
        )
      ) %>% 
      mutate(
        # Initial assignment for during follow-up (NA -> 0 first)
        !!bin_var := ifelse(
          is.na(.data[[bin_var]]),
          ifelse(
            elig_date >= .data[[start_var]] & elig_date <= .data[[end_var]],
            0, # assign 0 in current interval
            NA_real_
          ),
          .data[[bin_var]]
        )
      ) %>%
      arrange(.data[[start_var]]) %>%
      # Shift the 1 to next interval (lag)
      mutate(
        # create a temporary flag where event occurs
        event_here = ifelse(elig_date >= .data[[start_var]] & elig_date <= .data[[end_var]], 1, 0),
        # propagate 1 from next interval onwards
        !!bin_var := dplyr::lag(cummax(event_here), default = 0),
        # !!bin_var := ifelse(.data[[bin_var]] == 0 & event_cum == 1, 1, .data[[bin_var]]),
        # !!bin_var := zoo::na.locf(.data[[bin_var]], na.rm = FALSE), # forward fill
        # !!bin_var := replace_na(.data[[bin_var]], 0) # backward fill
      ) %>%
      ungroup() 
    # %>%
      # select(-elig_date, -event_here, -event_cum)
  })
  
  return(df_out)
}