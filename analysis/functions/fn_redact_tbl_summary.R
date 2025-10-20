################################################################################
# Custom made function to redact a tbl_summary object, it does the following:
# Take a tbl_summary() object as input and a threshold for redaction
# Redact categorical/dichotomous counts, including the total column total (second row N Total)
# Recalculate percentages for subgroup rows, based on the redacted column total
# Leave continuous variables untouched
# Return a clean tibble with renamed columns
## CAVE: inside uses the general midpoint rounding function: fn_roundmid_any
################################################################################
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

fn_redact_tbl_summary <- function(tbl_summary_obj, threshold = 6, # default threshold
                               col_names = c("stat_1", "stat_2"),
                               rename_cols = c("group1", "group2")) {
  
  # Step 0: extract table body as tibble
  tbl_df <- as_tibble(tbl_summary_obj$table_body)
  
  # Step 1: redact counts for categorical/dichotomous rows
  tbl_df <- tbl_df %>%
    mutate(
      across(
        all_of(col_names),
        ~ ifelse(
          var_type %in% c("categorical", "dichotomous"),
          map_chr(.x, function(cell) {
            if (is.na(cell) || cell == "") return(cell)
            txt <- str_replace_all(cell, ",", "")
            txt <- str_trim(txt)
            count <- suppressWarnings(as.numeric(str_extract(txt, "^[0-9]+")))
            if (is.na(count)) return(cell)
            fn_roundmid_any(count, threshold) %>% as.character()
          }),
          .x
        )
      )
    ) %>%
    # Step 2: extract numeric counts
    mutate(
      stat_1_count = suppressWarnings(as.numeric(!!sym(col_names[1]))),
      stat_2_count = suppressWarnings(as.numeric(!!sym(col_names[2])))
    )
  
  # Step 3: get redacted totals as scalars
  total_stat_1 <- tbl_df %>%
    filter(var_label == "Total N" & var_type == "categorical" & row_type == "level") %>%
    pull(stat_1_count) %>% as.numeric()
  
  total_stat_2 <- tbl_df %>%
    filter(var_label == "Total N" & var_type == "categorical" & row_type == "level") %>%
    pull(stat_2_count) %>% as.numeric()
  
  # Step 4: recalc percentages for non-total categorical/dichotomous rows
  tbl_df <- tbl_df %>%
    mutate(
      !!col_names[1] := ifelse(
        var_type %in% c("categorical", "dichotomous") & !var_label %in% c("Total N", "N") & !is.na(stat_1_count),
        paste0(formatC(stat_1_count, format = "f", digits = 0, big.mark = ","),
               " (", round(100 * stat_1_count / total_stat_1, 1), "%)"),
        !!sym(col_names[1])
      ),
      !!col_names[2] := ifelse(
        var_type %in% c("categorical", "dichotomous") & !var_label %in% c("Total N", "N") & !is.na(stat_2_count),
        paste0(formatC(stat_2_count, format = "f", digits = 0, big.mark = ","),
               " (", round(100 * stat_2_count / total_stat_2, 1), "%)"),
        !!sym(col_names[2])
      )
    ) %>%
    select(-stat_1_count, -stat_2_count) %>%
    # Step 5: optional rename for output clarity
    rename(!!rename_cols[1] := !!sym(col_names[1]),
           !!rename_cols[2] := !!sym(col_names[2]))
  
  return(tbl_df)
}
