####
# Custom-made function to extract confidence intervals from bootstrap object, for a specific timepoint, and keep a column with the original point estimates (without CI)
#### 
fn_extract_ci_boot <- function(boot_obj, index) {
  ci <- boot.ci(boot_obj, conf = 0.95, type = "perc", index = index)
  if (!is.null(ci$percent)) {
    return(c(ci$percent[4], ci$percent[5])) # Lower and Upper CI for the bootstrapped values
  } else {
    return(c(NA, NA)) # for the original value
  }
}