####
# Custom-made function to add dynamic time-updated covariates to person-interval data and adhere to specific rules:
#### 
# rule (a) If cov_date is in same interval as treatment info change but cov_date > treatment info change, then flag the one closest to end_date_month to be moved ("MOVER") to the next month (and discard other cov_date that are in same interval and with cov_date > treatment info change)
# rule (b) If cov_date is in same interval as treatment info change but cov_date <= treatment info change, then flag the one on or closest to treatment info change to keep ("KEEPER") in the same month (and discard the others that are also in same interval as treatment info change and cov_date <= treatment info change)
## (btw, if at the same time, there were also cov_date > treatment info change in same interval, they will have been dealt with in rule (a) already)
# rule (c) If cov_date is not in same interval as treatment info change, then flag the one closest to end_date_month to be moved ("MOVER") to the next month (and discard all others)
# Then, shift all "MOVERS" to the next month, and keep all "KEEPERS" in the original month

## Then, deal with edge cases: 
# (i) We may have months whereby we moved a "MOVER" covariate due to rule (c) or rule (a) into a month with a "KEEPER" covariate (rule b) => more than 1 covariate info in same person-interval
## => Solution: "KEEPER" beats "MOVER", i.e., the one that was flagged as "KEEPER" (because it contains time-update covariate info JUST before treatment change) is more important and wins
### (Side-note: In our case, we only have treatment update once, so a move due to rule (a) ending up in a month with a "KEEPER" covariate (rule b) is impossible - but the function/rules work with more complex dynamic treatment updates, too)
## (ii) We may have the scenario whereby treatment update info and covariate update info change on same date (and no time-stamp to differentiate)
## => Solution: Regard these covariate info as measured BEFORE treatment info change (see rule b above). Alternatively (sensitivity analyses), adapt rule a and b to regard them as measured AFTER treatment info change

fn_assign_dynamic_tu_cov <- function(data, 
                                     patient_id_col,
                                     variable_col, 
                                     date_col, 
                                     treat_date_col,
                                     month_col,
                                     start_date_col,
                                     end_date_col) {
  # Capture column names
  patient_id_col <- rlang::ensym(patient_id_col)
  variable_col <- rlang::ensym(variable_col)
  date_col <- rlang::ensym(date_col)
  treat_col <- rlang::ensym(treat_date_col)
  month_col <- rlang::ensym(month_col)
  start_col <- rlang::ensym(start_date_col)
  end_col <- rlang::ensym(end_date_col)
  
  data %>%
    # Group by patient, month, and variable to handle multiple measurements of the same type within one month per person
    group_by(!!patient_id_col, !!variable_col, !!month_col) %>%
    arrange(!!date_col, .by_group = TRUE) %>%
    
    # time_diff_to_end can be calculated irrespective if there is treatment change in that month, but 
    # time_diff_to_treat needs to be calculate ONLY for covariate date_col that happens before treat_col in same interval, see use further down
    # function handles treat_col = NA correctly: 
    ## No KEEPER is flagged for NA treatment months
    ## MOVERS are still considered based on rule (c)
    mutate(
      time_diff_to_end = as.numeric(!!end_col - !!date_col),
      time_diff_to_treat_before = if_else(!!date_col <= !!treat_col,
                                          as.numeric(!!date_col - !!treat_col),
                                          NA_real_)
    ) %>%
    
    mutate(
      # flag all (rows) that are in same interval as our treatment information change/update date variable
      is_treat_month = case_when(
        !is.na(!!treat_col) &
          !!treat_col >= !!start_col & !!treat_col <= !!end_col &
          !!date_col >= !!start_col & !!date_col <= !!end_col ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    
    # Precompute groupwise minima to avoid rowwise min() warnings
    mutate(
      min_time_diff_to_end = suppressWarnings(min(abs(time_diff_to_end), na.rm = TRUE)),
      min_time_diff_to_treat_before = suppressWarnings(min(abs(time_diff_to_treat_before), na.rm = TRUE))
    ) %>%
    
    mutate(
      # rule (a)
      is_treat_month_after_treat_move = case_when(
        is_treat_month & !!date_col > !!treat_col &
          abs(time_diff_to_end) == min_time_diff_to_end ~ TRUE,
        TRUE ~ FALSE
      ),
      # rule (b)
      is_treat_month_before_treat_keep = case_when(
        is_treat_month & !!date_col <= !!treat_col &
          abs(time_diff_to_treat_before) == min_time_diff_to_treat_before ~ TRUE,
        TRUE ~ FALSE
      ),
      # rule (c)
      is_not_treat_month_move = case_when(
        !is_treat_month &
          abs(time_diff_to_end) == min_time_diff_to_end ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    
    ungroup() %>%
    # discard all those that are not flagged as "MOVER" or "KEEPER"
    filter(is_not_treat_month_move | is_treat_month_before_treat_keep | is_treat_month_after_treat_move) %>%
    
    mutate(
      # shift the "MOVERS" to next month, keep the "KEEPERS" with the original month
      new_month = case_when(
        is_not_treat_month_move | is_treat_month_after_treat_move ~ (!!month_col) + 1,
        TRUE ~ (!!month_col)
      )
    ) %>%
    
    group_by(!!patient_id_col, !!variable_col, new_month) %>%
    # dealt with edge case (i): Due to the move, we now may have instances whereby we moved a "MOVER" one into a month with a "KEEPER" one.
    # Use "arrange", so that TRUE values for is_cens_month_before_cens_keep come first and then keep only those ones ("KEEPER" beats "MOVER")
    arrange(desc(is_treat_month_before_treat_keep)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # Drop helper columns
    select(-c(time_diff_to_end,
              time_diff_to_treat_before,
              min_time_diff_to_end,
              min_time_diff_to_treat_before,
              is_treat_month,
              is_treat_month_after_treat_move,
              is_treat_month_before_treat_keep,
              is_not_treat_month_move))
}
