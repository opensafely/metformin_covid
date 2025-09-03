####
# Custom-made function to add stable time-udpated covariates to person-interval data and adhere to specific rules:
#### 
# Rule a–c:
## If covariate date is NA -> all intervals = 0
## If covariate date < min follow-up start -> all intervals = 1
## If covariate date > max follow-up end -> all intervals = 0

# Rule d (during follow-up; covariate event not in same interval as treatment information change):
## Mark 0 up to and including the interval of the event,
## Then assign 1 from the next interval onwards (via lag(cummax())).

# Rule e (during follow-up; covariate event in same interval as treatment information change):
## If cov_date <= treat_date, then mark *current* interval as 1. Otherwise, keep as-is.

fn_assign_stable_tu_cov <- function(df, date_vars, start_var, end_var, cens_date_var) {
  
  df_out <- df
  
  walk(date_vars, function(var) {
    # create binary flagging variable, based on _date_ naming
    bin_var <- str_replace(var, "_date_", "_bin_")
    
    df_out <<- df_out %>%
      group_by(patient_id) %>%
      arrange(.data[[start_var]], .by_group = TRUE) %>%
      
      mutate(
        cov_date = .data[[var]],
        min_start = min(.data[[start_var]], na.rm = TRUE),
        max_end   = max(.data[[end_var]], na.rm = TRUE),
        
        # rules a–c
        !!bin_var := case_when(
          is.na(cov_date) ~ 0, 
          cov_date < min_start ~ 1, 
          cov_date > max_end   ~ 0, 
          TRUE ~ NA_real_   
        )
      ) %>%
      
    # rule d
    mutate(
      event_here = ifelse(!is.na(cov_date) &
                            cov_date >= .data[[start_var]] &
                            cov_date <= .data[[end_var]], 1, 0),
      
      cov_flag_temp = dplyr::lag(cummax(event_here), default = 0),
      
      !!bin_var := ifelse(is.na(.data[[bin_var]]), cov_flag_temp, .data[[bin_var]])
    ) %>%
      
    # rule e
    mutate(
      !!bin_var := case_when(
        # covariate in same interval as censoring
        !is.na(cov_date) & cov_date >= .data[[start_var]] & cov_date <= .data[[end_var]] &
          !is.na(.data[[cens_date_var]]) & .data[[cens_date_var]] >= .data[[start_var]] & .data[[cens_date_var]] <= .data[[end_var]] &
          cov_date <= .data[[cens_date_var]] ~ 1,  # if cov_date before or at censor date in same interval -> mark/overwrite current interval with 1
        TRUE ~ .data[[bin_var]]
      )
    ) %>%
      ungroup() %>%
      select(-cov_date, -min_start, -max_end, -event_here, -cov_flag_temp)
  })
  
  return(df_out)
}