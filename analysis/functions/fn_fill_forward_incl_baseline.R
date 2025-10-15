####
# Custom-made function to fill forward covariate information, based on:
# expanded interval dataset with time-updated covariate information in some intervals => fill forward taking the latest value into account
# If the first interval is missing, take the value from a separate corresponding baseline variable (e.g. the time-fixed cov_num_bmi_b to fill in the time-updated cov_num_bmi at baseline and then take it from there)
####

fn_fill_forward_incl_baseline <- function(x, baseline_col_name) {
  # pick() returns a data frame of the current group; select first row
  base_val <- pick(all_of(baseline_col_name))[[1]][1]  # <<-- take first row only
  
  first_non_na <- which(!is.na(x))[1]
  
  if (is.na(first_non_na)) {
    if (!is.na(base_val)) {
      rep(base_val, length(x))
    } else {
      x
    }
  } else {
    masked <- x
    if (!is.na(base_val)) {
      masked[1:(first_non_na - 1)] <- base_val
    } else {
      masked[1:(first_non_na - 1)] <- NA
    }
    
    # forward fill
    tidyr::fill(data.frame(masked), masked, .direction = "down")$masked
  }
}