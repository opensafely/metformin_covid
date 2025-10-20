####
# Custom-made function to calculate mean difference or/and standardized mean difference BEFORE any weighting
# In case both groups are NA (e.g. as in dummy data): skip
# If no levels exist or smd_levels is empty: set to NA.

fn_smd_before <- function(x) {
  t0 <- treat_b0[[x]]
  t1 <- treat_b1[[x]]
  
  # Handle binary/logical variables
  if (is.logical(t0) || is.logical(t1)) {
    t0 <- factor(t0, levels = c(FALSE, TRUE))  # Convert the logicals to factors
    t1 <- factor(t1, levels = c(FALSE, TRUE))  # Convert the logicals to factors
  }
  
  # Checks if either t0 or t1 contains ONLY missing values; if so, don't calculate SMD
  if (all(is.na(t0)) || all(is.na(t1))) {
    result <- c(var = x, type = "empty", smd = NA)
    
  } else if (is.numeric(t0)) {
    # Numeric: SMD
    pooled_sd <- sqrt((var(t0, na.rm = TRUE) + var(t1, na.rm = TRUE)) / 2)
    smd <- (mean(t1, na.rm = TRUE) - mean(t0, na.rm = TRUE)) / pooled_sd # divided by pooled standard deviation to get SMD
    result <- c(var = x, type = "numeric", smd = smd)
    
  } else if (is.factor(t0) || is.character(t0)) {
    # Categorical/Factor or Binary/Logical
    # Categorical: compute SMD for each level, report max! p_pool: the average of p0 and p1.
    levels_all <- union(levels(factor(t0)), levels(factor(t1)))
    smd_levels <- c()
    for (lvl in levels_all) {
      p0 <- mean(t0 == lvl, na.rm = TRUE)
      p1 <- mean(t1 == lvl, na.rm = TRUE)
      p_pool <- (p0 + p1) / 2
      # Avoid divide by zero (-> SMD set to 0)
      if (p_pool == 0 || p_pool == 1) {
        smd_lvl <- 0
      } else {
        smd_lvl <- (p1 - p0) / sqrt(p_pool * (1 - p_pool)) # otherwise use formula
      }
      smd_levels <- c(smd_levels, smd_lvl)
    }
    # any SMD calculated for each level of the categorical variable?
    if (length(smd_levels) == 0 || all(is.na(smd_levels))) {
      # If no valid SMD values exist (either because the vector is empty or all values are NA), then assign NA to max_smd...
      max_smd <- NA
    } else {
      # ...otherwise use max
      max_smd <- max(abs(smd_levels), na.rm = TRUE)
    }
    result <- c(var = x, type = "categorical", smd = max_smd)
    
  } else {
    result <- c(var = x, type = "unknown", smd = NA)
  }
  
  return(result)
}


