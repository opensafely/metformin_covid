####
# Custom-made function to fill forward covariate information, based on:
# expanded interval dataset with time-updated covariate information in some intervals => fill forward taking the latest value into account
# If the first interval is missing, take the value from a separate corresponding baseline variable (e.g. the time-fixed cov_num_bmi_b to fill in the time-updated cov_num_bmi at baseline and then take it from there)
####

fn_fill_forward_incl_baseline <- function(x, baseline_col_name) {
  # pick() returns a data frame of the current group (i.e. current covariate within a person) 
  base_val <- pick(all_of(baseline_col_name))[[1]][1]  # select first row, and convert to vector
  
  first_non_na <- which(!is.na(x))[1]
  
  # first deal with case where no observed values exist in x
  if (is.na(first_non_na)) {
    if (!is.na(base_val)) {
      rep(base_val, length(x)) # use baseline value (if baseline missing, then all NA)
    } else {
      x
    }
  } else {
    # in the case where at least one observed value exists in x
    masked <- x
    if (!is.na(base_val)) {
      masked[1:(first_non_na - 1)] <- base_val # fills all missing values before the first observed value with the baseline.
    } else {
      masked[1:(first_non_na - 1)] <- NA
    }
    tidyr::fill(data.frame(masked), masked, .direction = "down")$masked # forward fills missing values in the vector.
  }
}