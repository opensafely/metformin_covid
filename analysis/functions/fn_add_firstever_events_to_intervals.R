####
# Custom-made function to add one-time first-ever events to (monthy/weekly) interval person/month dataframes

# RULES:
# a) if event date is not NA and happened before min start date intervals, then assign 1 (in corresponding flag variable) to all person-intervals
# b) if event date is not NA and happened after max end date intervals, then assign 0 (in corresponding flag variable) to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 (in corresponding flag variable) to all person-intervals
# d) if event date happened during follow-up (i.e. existing interval date), then assign 1 (in corresponding flag variable) to corresponding person-interval, and 0 to all person-intervals before, and 1 to all person-intervals after

fn_add_firstever_events_to_intervals <- function(df, date_vars, start_var, end_var) {
  
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
          elig_date < min(.data[[start_var]], na.rm = TRUE) ~ 1,
          elig_date >= min(.data[[start_var]], na.rm = TRUE) & elig_date <= max(.data[[end_var]], na.rm = TRUE) ~ NA_real_,  # assign NA for now, deal with it in a second step
          elig_date > max(.data[[end_var]], na.rm = TRUE) ~ 0
        )
      ) %>% 
      mutate(
        !!bin_var := ifelse(
          is.na(.data[[bin_var]]),
          ifelse(
            elig_date >= .data[[start_var]] & elig_date <= .data[[end_var]],
            1,
            NA_real_
          ),
          .data[[bin_var]]
        )
      ) %>%
      arrange(.data[[start_var]]) %>%
      mutate(
        !!bin_var := zoo::na.locf(.data[[bin_var]], na.rm = FALSE), # assign 1 forward
        !!bin_var := replace_na(.data[[bin_var]], 0) # assign 0 backwards
      ) %>%
      ungroup() %>%
      select(-elig_date)
  })
  
  return(df_out)
}


