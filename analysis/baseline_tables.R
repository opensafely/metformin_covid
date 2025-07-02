####
## This script does the following:
# 1. Import released midpoint rounded disclosure-safe output data
# 2. Create descriptive baseline tables 
# 3. Save all output
####

# Import libraries and user functions -------------------------------------
library('tidyverse')
library('here')
library('gt')

# Create directories for output -------------------------------------------
fs::dir_create(here::here("output", "data_description"))

# Import the data ---------------------------------------------------------
print('Import the data')
df_main <- read_csv(here("output", "data_description", "table1_main_midpoint6.csv"))
df_death_ltfu1 <- read_csv(here("output", "data_description", "table1_death_ltfu1_midpoint6.csv"))
df_death_ltfu2 <- read_csv(here("output", "data_description", "table1_death_ltfu2_midpoint6.csv"))

# Create baseline table ---------------------------------------------------
print('Create baseline table')

# Step 1: Format the stats and create stats labels
print('Step 1: Format the stats and create stats labels')
stats_labels <- function(df) {
  df %>%
    mutate(var_label = case_when(is.na(var_label) ~ "Total N", TRUE ~ var_label)) %>%
    mutate(
      p = p * 100,
      p = round(p, 2),
      median = round(median, 2),
      p25 = round(p25, 2),
      p75 = round(p75, 2),
      mean = round(mean, 2),
      sd = round(sd, 2),
      stat_label = case_when(
        !is.na(median) ~ str_glue("{median} ({p25}, {p75}); {mean} ({sd})"),
        !is.na(n) ~ str_glue("{n} ({p}%)"),
        TRUE ~ ""
      )
    )
}
df_main <- stats_labels(df_main)
df_death_ltfu1 <- stats_labels(df_death_ltfu1)
df_death_ltfu2 <- stats_labels(df_death_ltfu2)

# Step 2: Create 'Unknown' level for each variable and each level
print('Step 2: Create Unknown level for each variable and each level')
add_unknown_levels <- function(df) {
  unknown_levels <- df %>%
    group_by(var_label, by) %>%
    summarise(
      variable_levels = "Unknown",
      stat_label = as.character(first(N_miss)),
      .groups = "drop"
    )
  
  df %>%
    bind_rows(unknown_levels) %>%
    arrange(var_label, variable_levels, by)
}
df_main <- add_unknown_levels(df_main)
df_death_ltfu1 <- add_unknown_levels(df_death_ltfu1)
df_death_ltfu2 <- add_unknown_levels(df_death_ltfu2)

# Step 3: Reshape the data to wide format
print('Step 3: Reshape the data to wide format')
reshape_to_wide <- function(df) {
  df %>% 
    select(var_label, variable_levels, by, stat_label) %>% 
    distinct() %>% # Ensure unique rows before pivoting
    pivot_wider(
      names_from = by, 
      values_from = stat_label,
      values_fn = list(stat_label = first) # Take the first value if duplicates exist
    )
}
df_main <- reshape_to_wide(df_main)
df_death_ltfu1 <- reshape_to_wide(df_death_ltfu1)
df_death_ltfu2 <- reshape_to_wide(df_death_ltfu2)

# Step 4: Order the variables and levels within variables
print('Step 4: Order the variables and levels within variables')
custom_order_main <- c("Total N", "Age", "Age groups", "Sex", "Ethnicity", "Deprivation", "Region", "Rural/urban", "Smoking status",
                  "Care/nursing home resident", "Healthcare worker", "Consultation rate in previous year", 
                  "Body Mass Index > 30 kg/m^2", "HbA1c in mmol/mol", "TC/Chol ratio", "HbA1c categories in mmol/mol", "TC/HDL ratio categories", 
                  
                  "History of acute myocardial infarct", "History of stroke", "History of other arterial embolism", "History of venous thromboembolism",
                  "History of heart failure", "History of angina pectoris", "History of dementia", "History of cancer", "History of arterial hypertension",
                  "History of depression", "History of COPD", "History of liver disease", "History of CKD", "History of PCOS", "History of prediabetes",
                  "Diabetes complication", 
                  
                  "Calendar period of T2DM diagnosis",

                  "COVID hosp or death", "COVID hosp", "COVID death", "Any covid diagnosis, pos test or hosp",

                  "Any Long COVID diagnosis",
                  "Any Long COVID or Viral Fatigue diagnosis",
                  "Any Viral Fatigue diagnosis",
                  "Any death after landmark",
                  "Any LTFU after landmark",
                  
                  "INT: Any metformin prescription within 6m prior to pandemic start", 
                  "CONT: Any metformin start in control"
                  )
custom_order_death_ltfu <- c("Total N", 
                              "Metformin treatment",
                              "Age", "Age groups", "Sex", "Ethnicity", "Deprivation", "Region", "Rural/urban", "Smoking status",
                       "Care/nursing home resident", "Healthcare worker", "Consultation rate in previous year", 
                       "Body Mass Index > 30 kg/m^2", "HbA1c in mmol/mol", "TC/Chol ratio", "HbA1c categories in mmol/mol", "TC/HDL ratio categories", 
                       
                       "History of acute myocardial infarct", "History of stroke", "History of other arterial embolism", "History of venous thromboembolism",
                       "History of heart failure", "History of angina pectoris", "History of dementia", "History of cancer", "History of arterial hypertension",
                       "History of depression", "History of COPD", "History of liver disease", "History of CKD", "History of PCOS", "History of prediabetes",
                       "Diabetes complication", 
                       
                       "Calendar period of T2DM diagnosis",
                       
                       "COVID hosp or death", "COVID hosp", "COVID death", "Any covid diagnosis, pos test or hosp",
                       
                       "Any Long COVID diagnosis",
                       "Any Long COVID or Viral Fatigue diagnosis",
                       "Any Viral Fatigue diagnosis",
                       "Any death after landmark",
                       "Any LTFU after landmark",
                       
                       "INT: Any metformin prescription within 6m prior to pandemic start", 
                       "CONT: Any metformin start in control"
)

reorder_baseline_table <- function(df, custom_order) {
  df %>%
    mutate(
      # Ensure var_label is ordered with custom_order first, others after
      var_label = factor(
        var_label,
        levels = c(custom_order, setdiff(unique(var_label), custom_order)),
        ordered = TRUE
      )
    ) %>%
    arrange(var_label) %>% # Unspecified variables are kept but pushed to the end
    
    # Handle NA values in variable_levels (temporarily replace with "NA")
    mutate(variable_levels = ifelse(is.na(variable_levels), "NA", variable_levels)) %>%
    
    # Move 'Unknown' to the bottom of 'variable_levels' within each 'var_label'
    mutate(
      variable_levels = factor(
        variable_levels,
        levels = c(setdiff(unique(variable_levels), "Unknown"), "Unknown")
      )
    ) %>%
    # Now, sort
    arrange(var_label, variable_levels) %>%
    ungroup()
}
df_main <- reorder_baseline_table(df_main, custom_order_main)
df_death_ltfu1 <- reorder_baseline_table(df_death_ltfu1, custom_order_death_ltfu)
df_death_ltfu2 <- reorder_baseline_table(df_death_ltfu2, custom_order_death_ltfu)

# Step 5: Drop "Unknown" rows where var_label has already Unknown as a level or not needed
print('Step 5: Drop "Unknown" rows where var_label has already Unknown as a level or not needed')
drop_unknown <- function(df) {
  df %>%
    filter(!(var_label == "Total N" & variable_levels == "Unknown")) %>%
    filter(!(grepl("History of", var_label) & variable_levels == "Unknown")) %>% 
    filter(!(grepl("metformin", var_label) & variable_levels == "Unknown")) %>%
    filter(!(grepl("COVID", var_label) & variable_levels == "Unknown")) %>%
    filter(!(grepl("Any", var_label) & variable_levels == "Unknown")) %>% 
    filter(!(var_label == "Diabetes complication" & variable_levels == "Unknown")) %>% 
    filter(!(var_label == "Calendar period of T2DM diagnosis" & variable_levels == "Unknown"))
  # filter(!(var_label == "Deprivation" & variable_levels == "Unknown")) %>%
  # filter(!(var_label == "Care/nursing home resident" & variable_levels == "Unknown")) %>% 
  # filter(!(var_label == "Healthcare worker" & variable_levels == "Unknown")) %>% 
  # filter(!(var_label == "Body Mass Index > 30 kg/m^2" & variable_levels == "Unknown"))
}
df_main <- drop_unknown(df_main)
df_death_ltfu1 <- drop_unknown(df_death_ltfu1)
df_death_ltfu2 <- drop_unknown(df_death_ltfu2)
  
# Step 6: Only keep the variable label for the first row of each variable group
print('Step 6: Only keep the variable label for the first row of each variable group')
first_row_label <- function(df) {
  df %>%
    mutate(
      var_label = if_else(
        row_number() == 1 | var_label != lag(var_label),
        var_label,
        factor(NA, levels = levels(var_label), ordered = is.ordered(var_label))
      )
    )
}
df_main <- first_row_label(df_main)
df_death_ltfu1 <- first_row_label(df_death_ltfu1)
df_death_ltfu2 <- first_row_label(df_death_ltfu2)

# Step 7: Hide NAs where not needed
print('Step 7: Hide NAs where not needed')
hide_NA <- function(df) {
  df %>%
    mutate(variable_levels = recode(variable_levels, "NA" = "")) %>% 
    mutate(var_label = replace_na(as.character(var_label), ""))
}
df_main <- hide_NA(df_main)
df_death_ltfu1 <- hide_NA(df_death_ltfu1)
df_death_ltfu2 <- hide_NA(df_death_ltfu2)

# Step 8: Save the underlying data as csv
print('Step 8: Save the underlying data as csv')
tbl_csv_main <- df_main %>%
  select(var_label, variable_levels, `Metformin mono`, Nothing)
tbl_csv_death_ltfu1 <- df_death_ltfu1 %>%
  select(var_label, variable_levels, `Alive at landmark`, `Died until landmark`)
tbl_csv_death_ltfu2 <- df_death_ltfu2 %>%
  select(var_label, variable_levels, `Alive and in care at pandemic start`, `Died or LTFU between landmark and pandemic start`)

# Step 9: Create the gt table
print('Step 9: Create the gt table')
tbl_gt_main <- df_main %>%
  select(var_label, variable_levels, `Metformin mono`, Nothing) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    var_label = "Characteristic",
    variable_levels = "",  # Empty label
    `Metformin mono` = "Metformin mono",
    Nothing = "Nothing"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )
tbl_gt_death_ltfu1 <- df_death_ltfu1 %>%
  select(var_label, variable_levels, `Alive at landmark`, `Died until landmark`) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    var_label = "Characteristic",
    variable_levels = "",  # Empty label
    `Alive at landmark` = "Alive at landmark",
    `Died until landmark` = "Died until landmark"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )
tbl_gt_death_ltfu2 <- df_death_ltfu2 %>%
  select(var_label, variable_levels, `Alive and in care at pandemic start`, `Died or LTFU between landmark and pandemic start`) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    var_label = "Characteristic",
    variable_levels = "",  # Empty label
    `Alive and in care at pandemic start` = "Alive and in care at pandemic start",
    `Died or LTFU between landmark and pandemic start` = "Died or LTFU between landmark and pandemic start"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )

# Save baseline tables ---------------------------------------------------
print('Save baseline tables')
write.csv(tbl_csv_main, file = here::here("output", "data_description", "tbl_csv_main.csv"))
write.csv(tbl_csv_death_ltfu1, file = here::here("output", "data_description", "tbl_csv_death_ltfu1.csv"))
write.csv(tbl_csv_death_ltfu2, file = here::here("output", "data_description", "tbl_csv_death_ltfu2.csv"))
