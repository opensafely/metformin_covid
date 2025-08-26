####
# Custom-made function to add dynamic time-updated covariates and adhere to specific rules:
#### 
# Now, if there are several of the same dynamic time-updated covariate events happening in the same interval/month (e.g. HbA1c several times recorded/measured in same month):
## => carefully choose the relevant one, according to the following rules, taking our treatment information change/update variable (cens_date_metfin_start_cont) into account: 
### a) If cov_date is in same interval as cens_date_metfin_start_cont and cov_date > cens_date_metfin_start_cont, then flag the one closest to end_date_month to be moved to the next month (and discard the others that are in same interval as cens_date_metfin_start_cont & with cov_date > cens_date_metfin_start_cont)
### b) If cov_date is in same interval as cens_date_metfin_start_cont and cov_date <= cens_date_metfin_start_cont, then flag the one on or closest to cens_date_metfin_start_cont to keep (and discard the others that are also in same interval as cens_date_metfin_start_cont & cov_date <= cens_date_metfin_start_cont)
#### If at the same time, there were also cov_date > cens_date_metfin_start_cont in same interval, they will have been dealt with in rule (a) already
### c) If cov_date is not in same interval as cens_date_metfin_start_cont, then flag the one closest to end_date_month to be moved to the next month (and discard all others)
### d) Then, shift all covariate info flagged as "move" to the next month, while keeping all flagged as "keep" in the original month

##### Then, deal with edge cases: 
### (i) We may have months whereby we moved a covariate due to rule (c) or rule (a) into a month with a "keep" covariate (rule b) => >1 covariate info update in same person-interval
#### => Solution: "keep" beats "move", i.e., the one that was flagged as "keep" (because it contains time-update covariate info JUST before treatment change) is more important and wins
##### Side-note: In our case, we only have treatment update once, so a move due to rule (a) ending up in a month with a "keep" covariate (rule b) is not possible - but the function/rules work with more complex dynamic treatment updates, too!
### (ii) We may have the scenario whereby treatment update info and covariate update info change on same date (and no time-stamp to differentiate)
#### => Solution: Regard these covariate info as measured BEFORE treatment info change (see rule b above). Alternatively (sensitivity analyses), simply adapt rule a and b to regard them as measured AFTER treatment info change


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
    
    mutate(
      time_diff_to_end = as.numeric(!!end_col - !!date_col),
      time_diff_to_treat = as.numeric(!!date_col - !!treat_col)
    ) %>%
    
    mutate(
      # flag all (rows) that are in same interval as our treatment information change/update variable (cens_date_metfin_start_cont)
      is_treat_month = case_when(
        !is.na(!!treat_col) &
          !!treat_col >= !!start_col & !!treat_col <= !!end_col &
          !!date_col >= !!start_col & !!date_col <= !!end_col ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    
    # Precompute groupwise minima to avoid rowwise min() warnings
    mutate(
      min_time_diff_to_end = suppressWarnings(min(abs(time_diff_to_end),   na.rm = TRUE)),
      min_time_diff_to_treat = suppressWarnings(min(abs(time_diff_to_treat), na.rm = TRUE))
    ) %>%
    
    mutate(
      # rule a
      is_treat_month_after_treat_move = case_when(
        is_treat_month & !!date_col > !!treat_col &
          abs(time_diff_to_end) == min_time_diff_to_end ~ TRUE,
        TRUE ~ FALSE
      ),
      # rule b
      is_treat_month_before_treat_keep = case_when(
        is_treat_month & !!date_col <= !!treat_col &
          abs(time_diff_to_treat) == min_time_diff_to_treat ~ TRUE,
        TRUE ~ FALSE
      ),
      # rule c
      is_not_treat_month_move = case_when(
        !is_treat_month &
          abs(time_diff_to_end) == min_time_diff_to_end ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    
    ungroup() %>%
    # discard all those that are not flagged as "move" or "keep"
    filter(is_not_treat_month_move | is_treat_month_before_treat_keep | is_treat_month_after_treat_move) %>%
    
    mutate(
      # shift the "move" ones to next month, keep the "keep" ones with the initial month
      new_month = case_when(
        is_not_treat_month_move | is_treat_month_after_treat_move ~ (!!month_col) + 1,
        TRUE ~ (!!month_col)
      )
    ) %>%
    
    group_by(!!patient_id_col, !!variable_col, new_month) %>%
    # dealt with edge case (i): Due to the move, we now may have instances whereby we moved a "move" one into a month with a "keep" one.
    # Use "arrange", so that TRUE values for is_cens_month_before_cens_keep come first and then keep only those ones ("keep" beats "move")
    arrange(desc(is_treat_month_before_treat_keep)) %>%
    slice(1) %>%
    ungroup() %>%
    
    # Drop helper columns
    select(-c(time_diff_to_end,
              time_diff_to_treat,
              min_time_diff_to_end,
              min_time_diff_to_treat,
              is_treat_month,
              is_treat_month_after_treat_move,
              is_treat_month_before_treat_keep,
              is_not_treat_month_move))
}
