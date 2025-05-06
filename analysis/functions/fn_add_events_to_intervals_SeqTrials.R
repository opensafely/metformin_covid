####
# Custom-made function to add events to (monthy/weekly) interval person/month dataframes

# RULES:
# a) if event date is recorded and happening before min start date intervals, then assign 1 to corresponding bin variable to all person-intervals
# b) if event date is recorded and happening after max end date intervals, then assign 0 to corresponding bin variable to all person-intervals
# c) if no event date is recorded (date variable == NA), then assign 0 to corresponding bin variable to all person-intervals
# d) if event date is recorded during follow-up, then assign 1 to corresponding bin variable in corresponding interval and 0 to all intervals before and 1 to all intervals after

fn_add_events_to_intervals <- function(df, date_vars, id_var, start_var, end_var) {
  
  df_out <- df
  
  walk(date_vars, function(var) {
    # Dynamically create binary variable name by replacing the _date_ with a corresponding _bin_
    bin_var <- str_replace(var, "_date_", "_bin_")
    
    df_out <<- df_out %>%
      group_by(.data[[id_var]]) %>%
      mutate(
        elig_date = first(na.omit(.data[[var]])),
        !!bin_var := case_when(
          is.na(elig_date) ~ 0,
          elig_date < min(.data[[start_var]]) ~ 1,
          elig_date >= min(.data[[start_var]]) & elig_date <= max(.data[[end_var]]) ~ NA_real_, # assign NA for now, deal with it in a second step
          elig_date > max(.data[[end_var]]) ~ 0
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


