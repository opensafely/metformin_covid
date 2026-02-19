####
# Custom-made function to truncate
#### 

# Args:
#   df: The data frame containing the column to be truncated.
#   col_name: A string specifying the name of the column to truncate.
#   lower_p: The lower percentile value (e.g., 0.01 for the 1st percentile).
#   upper_p: The upper percentile value (e.g., 0.99 for the 99th percentile).
#
# Returns:
#   The original data frame with a new column that has the original
#   column name plus a "_trunc" suffix.

fn_truncate_by_percentile <- function(df, col_name, lower_p = 0.01, upper_p = 0.99) {
  # Check if the column exists in the data frame
  if (!col_name %in% names(df)) {
    stop(paste("Column", col_name, "not found in the data frame."))
  }
  
  # Calculate the lower and upper percentile thresholds
  lower_threshold <- quantile(df[[col_name]], lower_p, na.rm = TRUE)
  upper_threshold <- quantile(df[[col_name]], upper_p, na.rm = TRUE)
  
  # Create a new column name with the "_trunc" suffix
  truncated_col_name <- paste0(col_name, "_trunc")
  
  # Create a new column for the truncated values
  df_truncated <- df
  df_truncated[[truncated_col_name]] <- df_truncated[[col_name]]
  
  # Truncate values below the lower threshold
  df_truncated[[truncated_col_name]][df_truncated[[truncated_col_name]] < lower_threshold] <- lower_threshold
  
  # Truncate values above the upper threshold
  df_truncated[[truncated_col_name]][df_truncated[[truncated_col_name]] > upper_threshold] <- upper_threshold
  
  # Return the data frame with the new truncated column
  return(df_truncated)
}