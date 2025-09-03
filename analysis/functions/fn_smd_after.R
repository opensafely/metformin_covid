####
# Custom-made function to calculate mean difference or/and standardized mean difference AFTER weighting
# In case both groups are NA (e.g. as in dummy data): skip
# If no levels exist or smd_levels is empty: set to NA.

fn_smd_after <- function(x, w_treat_0, w_treat_1) {
  # Subsetting the data based on the variable name (x)
  t0 <- treat_b0[[x]]
  t1 <- treat_b1[[x]]
  w0 <- w_treat_0  # Weights for the control group
  w1 <- w_treat_1  # Weights for the treated group
  
  # Handle binary/logical variables
  if (is.logical(t0) || is.logical(t1)) {
    t0 <- factor(t0, levels = c(FALSE, TRUE))  # Convert the logicals to factors
    t1 <- factor(t1, levels = c(FALSE, TRUE)) 
  }
  
  # Handle missing data
  if (all(is.na(t0)) || all(is.na(t1))) {
    result <- c(var = x, type = "empty", smd = NA)
  } else if (is.numeric(t0)) {
    # Weighted Numeric: SMD
    weighted_mean_t0 <- weighted.mean(t0, w0, na.rm = TRUE)
    weighted_mean_t1 <- weighted.mean(t1, w1, na.rm = TRUE)
    weighted_var_t0 <- sum(w0 * (t0 - weighted_mean_t0)^2, na.rm = TRUE) / sum(w0, na.rm = TRUE)
    weighted_var_t1 <- sum(w1 * (t1 - weighted_mean_t1)^2, na.rm = TRUE) / sum(w1, na.rm = TRUE)
    pooled_sd <- sqrt((weighted_var_t0 + weighted_var_t1) / 2)
    smd <- (weighted_mean_t1 - weighted_mean_t0) / pooled_sd
    result <- c(var = x, type = "numeric", smd = smd)
    
  } else if (is.factor(t0) || is.character(t0)) {
    # Weighted Categorical: compute SMD per level, report max
    levels_all <- union(levels(factor(t0)), levels(factor(t1)))
    smd_levels <- c()
    for (lvl in levels_all) {
      p0 <- weighted.mean(t0 == lvl, w0, na.rm = TRUE)
      p1 <- weighted.mean(t1 == lvl, w1, na.rm = TRUE)
      p_pool <- (p0 + p1) / 2
      # Avoid divide by zero
      if (p_pool == 0 || p_pool == 1) {
        smd_lvl <- 0
      } else {
        smd_lvl <- (p1 - p0) / sqrt(p_pool * (1 - p_pool))
      }
      smd_levels <- c(smd_levels, smd_lvl)
    }
    if (length(smd_levels) == 0 || all(is.na(smd_levels))) {
      max_smd <- NA
    } else {
      max_smd <- max(abs(smd_levels), na.rm = TRUE)
    }
    result <- c(var = x, type = "categorical", smd = max_smd)
    
  } else {
    result <- c(var = x, type = "unknown", smd = NA)
  }
  
  return(result)
}




