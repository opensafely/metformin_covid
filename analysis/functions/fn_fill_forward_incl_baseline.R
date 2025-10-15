####
# Custom-made function to fill forward covariate information, based on:
# expanded interval dataset with time-updated covariate information in some intervals => fill forward taking the latest value into account
# If the first interval is missing, take the value from a separate corresponding baseline variable (e.g. the time-fixed cov_num_bmi_b to fill in the time-updated cov_num_bmi at baseline and then take it from there)
####

fn_fill_forward_incl_baseline <- function(x, baseline_vec) {
  # x: time-updated vector | baseline_vec: corresponding baseline vector (same length as x, usually repeated per patient)
  
  base_val <- baseline_vec[1]
  
  # if first element is NA, fill it with baseline (if available)
  if (is.na(x[1]) && !is.na(base_val)) {
    x[1] <- base_val
  }
  
  # forward-fill | na.rm = FALSE keeps leading NAs if no baseline
  zoo::na.locf(x, na.rm = FALSE)
}