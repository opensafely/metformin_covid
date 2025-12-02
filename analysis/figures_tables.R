####
## This script does the following:
# 1. Import released midpoint rounded disclosure-safe output data
# 2. Create various tables and figures 
####

# Import libraries and user functions -----------------------------
library('dplyr')
library('tidyr')
library('readr')
# library('stringr')
library('here')
# library('DiagrammeR')
library('ggplot2')
# library('RColorBrewer')
library('gt')
library(forestplot)


# Import the data -------------------------------------------------
df_elig <- read_csv(here("output", "data_release_20250821", "data_description", "n_elig_excluded_midpoint6.csv"))
# df_qa <- read_csv(here("output", "data_release_20250303", "n_qa_excluded_midpoint6.csv"))
df_comp <- read_csv(here("output", "data_release_20250821", "data_description", "n_completeness_excluded_midpoint6.csv"))
df_elig_landmark <- read_csv(here("output", "data_release_20250821", "data_description", "n_restricted_midpoint6.csv"))
df_main <- read_csv(here("output", "data_release_20250821", "data_description", "table1_main_midpoint6.csv"))
df_ps <- read_csv(here("output", "data_release_20250821", "ps", "density_plot_untrimmed.csv"))
df_results_cox <- read_csv(here("output", "data_release_20250821", "make_output", "results_cox_midpoint6.csv"))


# Quality assurance flowchart -------------------------------------

# # Create node labels
# node_labels <- paste(df_qa$Variable, '\n', df_qa$Value)
# 
# # Create flowchart
# graph_qa <- grViz(
#   sprintf(
#     "digraph flowchart {
#       graph [layout = dot, rankdir = TB]
#       node [shape = box, style = filled, fontname = Arial, fontsize = 14]
# 
#       A [label = '%s', fillcolor = lightskyblue]
#       B [label = 'Exclusions\n%s\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       C [label = '%s', fillcolor = lightskyblue]
# 
#       A -> B
#       B -> C
#     }",
#     node_labels[1], node_labels[2],
#     node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
#     node_labels[8], node_labels[9]
#   )
# )

# Completeness criteria flowchart -----------------------------------

# # Create node labels
# node_labels <- paste(df_comp$Variable, '\n', df_comp$Value)
# 
# # Create flowchart
# graph_completeness <- grViz(
#   sprintf(
#     "digraph flowchart {
#       graph [layout = dot, rankdir = TB]
#       node [shape = box, style = filled, fontname = Arial, fontsize = 14]
# 
#       A [label = '%s', fillcolor = lightskyblue]
#       B [label = 'Exclusions\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       C [label = '%s', fillcolor = lightskyblue]
# 
#       A -> B
#       B -> C
#     }",
#     node_labels[1], node_labels[2],
#     node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
#     node_labels[8], node_labels[9]
#   )
# )


# Eligibility flowchart -----------------------------------

# # Create node labels
# node_labels <- paste(df_elig$Variable, '\n', df_elig$Value)
# 
# # Create flowchart
# graph_elig <- grViz(
#   sprintf(
#     "digraph flowchart {
#       graph [layout = dot, rankdir = TB]
#       node [shape = box, style = filled, fontname = Arial, fontsize = 14]
# 
#       A [label = '%s', fillcolor = lightskyblue]
#       B [label = '%s', fillcolor = lightyellow]
#       C [label = 'Exclusions 1\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       D [label = 'Exclusions 2: Events between baseline and landmark\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
#       E [label = '%s', fillcolor = lightskyblue]
# 
#       A -> B
#       B -> C
#       C -> D
#       D -> E
#     }",
#     node_labels[1], node_labels[2], 
#     node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
#     node_labels[8], node_labels[9], node_labels[10], node_labels[11], node_labels[12], 
#     node_labels[13], node_labels[14], node_labels[15], 
#     node_labels[16], node_labels[17], node_labels[18], node_labels[19], node_labels[20], node_labels[21],
#     node_labels[22]
#   )
# )


# Baseline table, based on output from table1.R ----------------------

# Step 1: Format the stats and create stats labels
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
# df_death_ltfu1 <- stats_labels(df_death_ltfu1)
# df_death_ltfu2 <- stats_labels(df_death_ltfu2)

# Step 2: Create 'Unknown' level for each variable and each level
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
# df_death_ltfu1 <- add_unknown_levels(df_death_ltfu1)
# df_death_ltfu2 <- add_unknown_levels(df_death_ltfu2)

# Step 3: Reshape the data to wide format
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
# df_death_ltfu1 <- reshape_to_wide(df_death_ltfu1)
# df_death_ltfu2 <- reshape_to_wide(df_death_ltfu2)

# Step 4: Order the variables and levels within variables
custom_order_main <- c("Total N", "Age", "Age groups", "Sex", "Ethnicity", "Deprivation", "Region", "Rural/urban", "Smoking status",
                       "Healthcare worker", "Consultation rate in previous year", 
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
                             "Healthcare worker", "Consultation rate in previous year", 
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
# df_death_ltfu1 <- reorder_baseline_table(df_death_ltfu1, custom_order_death_ltfu)
# df_death_ltfu2 <- reorder_baseline_table(df_death_ltfu2, custom_order_death_ltfu)

# Step 4.5: Reduce main table to main variables
df_main <- df_main %>% 
  filter(!grepl("INT: Any", var_label)) %>% 
  filter(!grepl("Any LTFU", var_label)) %>% 
  filter(!grepl("Any death", var_label)) %>% 
  filter(!grepl("Any Viral Fatigue", var_label)) %>% 
  filter(!grepl("Any Long COVID diagnosis", var_label))


# Step 5: Drop "Unknown" rows where var_label has already Unknown as a level or not needed
drop_unknown <- function(df) {
  df %>%
    filter(!(var_label == "Total N" & variable_levels == "Unknown")) %>%
    filter(!(grepl("History of", var_label) & variable_levels == "Unknown")) %>% 
    filter(!(grepl("metformin", var_label) & variable_levels == "Unknown")) %>%
    filter(!(grepl("COVID", var_label) & variable_levels == "Unknown")) %>%
    filter(!(grepl("Any", var_label) & variable_levels == "Unknown")) %>% 
    filter(!(var_label == "Diabetes complication" & variable_levels == "Unknown")) %>% 
    filter(!(var_label == "Calendar period of T2DM diagnosis" & variable_levels == "Unknown")) %>% 
    filter(!(var_label == "Healthcare worker" & variable_levels == "Unknown"))
  # filter(!(var_label == "Deprivation" & variable_levels == "Unknown")) %>%
  # filter(!(var_label == "Body Mass Index > 30 kg/m^2" & variable_levels == "Unknown"))
}
df_main <- drop_unknown(df_main)
# df_death_ltfu1 <- drop_unknown(df_death_ltfu1)
# df_death_ltfu2 <- drop_unknown(df_death_ltfu2)

# Step 6: Only keep the variable label for the first row of each variable group
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
# df_death_ltfu1 <- first_row_label(df_death_ltfu1)
# df_death_ltfu2 <- first_row_label(df_death_ltfu2)

# Step 7: Hide NAs where not needed
hide_NA <- function(df) {
  df %>%
    mutate(variable_levels = recode(variable_levels, "NA" = "")) %>% 
    mutate(var_label = replace_na(as.character(var_label), ""))
}
df_main <- hide_NA(df_main)
# df_death_ltfu1 <- hide_NA(df_death_ltfu1)
# df_death_ltfu2 <- hide_NA(df_death_ltfu2)

# Step 8: Save the underlying data as csv
tbl_csv_main <- df_main %>%
  select(var_label, variable_levels, `Metformin mono`, Nothing)
# tbl_csv_death_ltfu1 <- df_death_ltfu1 %>%
#   select(var_label, variable_levels, `Alive and in care at landmark`, `Died or LTFU until landmark`)
# tbl_csv_death_ltfu2 <- df_death_ltfu2 %>%
#   select(var_label, variable_levels, `Alive and in care at pandemic start`, `Died or LTFU between landmark and pandemic start`)

# Step 9: Create the gt table
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

print(tbl_gt_main)
# Save as RTF (Word-readable)
gtsave(tbl_gt_main, "baseline_characteristics.rtf")


# tbl_gt_death_ltfu1 <- df_death_ltfu1 %>%
#   select(var_label, variable_levels, `Alive and in care at landmark`, `Died or LTFU until landmark`) %>%
#   gt() %>%
#   tab_header(title = "Baseline Characteristics") %>%
#   cols_label(
#     var_label = "Characteristic",
#     variable_levels = "",  # Empty label
#     `Alive and in care at landmark` = "Alive and in care at landmark",
#     `Died or LTFU until landmark` = "Died or LTFU until landmark"
#   ) %>%
#   tab_options(
#     table.font.size = px(12)
#   )
# tbl_gt_death_ltfu2 <- df_death_ltfu2 %>%
#   select(var_label, variable_levels, `Alive and in care at pandemic start`, `Died or LTFU between landmark and pandemic start`) %>%
#   gt() %>%
#   tab_header(title = "Baseline Characteristics") %>%
#   cols_label(
#     var_label = "Characteristic",
#     variable_levels = "",  # Empty label
#     `Alive and in care at pandemic start` = "Alive and in care at pandemic start",
#     `Died or LTFU between landmark and pandemic start` = "Died or LTFU between landmark and pandemic start"
#   ) %>%
#   tab_options(
#     table.font.size = px(12)
#   )


# PS figure ----------------------
density_plot_untrimmed <- ggplot(df_ps, aes(x = dens_x, y = dens_y, color = group)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1), 
    minor_breaks = seq(0, 1, by = 0.05), 
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Propensity Score Density Plot; Untrimmed",
    x = "Propensity Score",
    y = "Density",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10)  # optional: make labels clearer
  )


# Results figures ----------------------
## Main results ----
df_results_cox <- df_results_cox %>% 
  filter(term == "days0_730")

df_results_main <- df_results_cox %>% 
  filter(!grepl("sg_", name)) %>% 
  filter(!grepl("_landmark", name)) %>% 
  filter(grepl("mdl_max_adj", model))

df_results_main$name <- recode(df_results_main$name,
                               "primary" = "Primary: COVID-19-related hospitalization or death",
                               "covid_event" = "Secondary: Any COVID-19-related event",
                               "longvirfat" = "Secondary: Long COVID",
                               "neg_control_pandemic" = "Negative control outcome (bone fracture)",
                               "pos_control_pandemic" = "Positive control outcome (diabetes-related death)",
)
desired_order <- c(
  "Primary: COVID-19-related hospitalization or death",
  "Secondary: Any COVID-19-related event",
  "Secondary: Long COVID",
  "Negative control outcome (bone fracture)",
  "Positive control outcome (diabetes-related death)"
)
df_results_main <- df_results_main %>%
  mutate(
    name = factor(name, levels = desired_order)
  ) %>%
  arrange(name)

# build forestplot
base_data <- tibble(
  mean = df_results_main$hr,
  lower = df_results_main$conf_low,
  upper = df_results_main$conf_high,
  name = as.character(df_results_main$name),
  # events_i = as.character(df_results_main$sarilumab),
  # events_c = as.character(df_results_main$tocilizumab),
  estimates = paste0(formatC(df_results_main$hr, format = "f", digits = 2), 
                     " (", formatC(df_results_main$conf_low, format = "f", digits = 2), 
                     " - ", formatC(df_results_main$conf_high, format = "f", digits = 2), ")"))
header <- tibble(
  name = "Outcome",
  # events_i = "Events Sarilumab\n n (%)",
  # events_c = "Events Tocilizumab\n n (%)",
  estimates = "aHR (95% CI)",
  mean = NA, lower = NA, upper = NA)

fp <- bind_rows(header, base_data)
font <- "sans"

# Create the forest plot
pdf(file = "Fp_results_main.pdf",
    width = 14,
    height = 7)

fp %>%
  forestplot(
    labeltext = c(name, 
                  # events_i, events_c, 
                  estimates),
    mean = mean, lower = lower, upper = upper,  # Numeric columns for the plot
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = font, cex = 1),
      ticks = gpar(cex = 0.88),
      xlab = gpar(cex = 0.88)
    ),
    graph.pos = 2, 
    hrzl_lines = list("2" = gpar(lty = 2),
                      "3" = gpar(lty = 2),
                      "4" = gpar(lty = 2),
                      "5" = gpar(lty = 2),
                      "6" = gpar(lty = 2),
                      "7" = gpar(lty = 2)),
    xlog = TRUE,
    # xticks = log(c(0.80, 1, 1.25)),
    lty.ci = 1,
    col = fpColors(
      box = "maroon4",
      line = "maroon1",
      summary = "magenta4",
      hrz_lines = "gray63"
    ),
    vertices = TRUE,
    xlab = NULL,
    zero = 1
  )

grid.text(
  "Higher outcome risk in no-treatment group < > Higher outcome risk in treatment group",
  x = unit(0.488, "npc"),
  y = unit(0.015, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()

## Minimally adjusted results ----
df_results_min <- df_results_cox %>% 
  filter(!grepl("sg_", name)) %>% 
  filter(!grepl("_landmark", name)) %>% 
  filter(!grepl("mdl_max_adj", model))

df_results_min$name <- recode(df_results_min$name,
                               "primary" = "Primary: COVID-19-related hospitalization or death",
                               "covid_event" = "Secondary: Any COVID-19-related event",
                               "longvirfat" = "Secondary: Long COVID",
                               "neg_control_pandemic" = "Negative control outcome (bone fracture)",
                               "pos_control_pandemic" = "Positive control outcome (diabetes-related death)",
)
desired_order <- c(
  "Primary: COVID-19-related hospitalization or death",
  "Secondary: Any COVID-19-related event",
  "Secondary: Long COVID",
  "Negative control outcome (bone fracture)",
  "Positive control outcome (diabetes-related death)"
)
df_results_min <- df_results_min %>%
  mutate(
    name = factor(name, levels = desired_order)
  ) %>%
  arrange(name)

# build forestplot
base_data <- tibble(
  mean = df_results_min$hr,
  lower = df_results_min$conf_low,
  upper = df_results_min$conf_high,
  name = as.character(df_results_min$name),
  # events_i = as.character(df_results_min$sarilumab),
  # events_c = as.character(df_results_min$tocilizumab),
  estimates = paste0(formatC(df_results_min$hr, format = "f", digits = 2), 
                     " (", formatC(df_results_min$conf_low, format = "f", digits = 2), 
                     " - ", formatC(df_results_min$conf_high, format = "f", digits = 2), ")"))
header <- tibble(
  name = "Outcome",
  # events_i = "Events Sarilumab\n n (%)",
  # events_c = "Events Tocilizumab\n n (%)",
  estimates = "aHR (95% CI)",
  mean = NA, lower = NA, upper = NA)

fp <- bind_rows(header, base_data)
font <- "sans"

# Create the forest plot
pdf(file = "Fp_results_min.pdf",
    width = 14,
    height = 7)

fp %>%
  forestplot(
    labeltext = c(name, 
                  # events_i, events_c, 
                  estimates),
    mean = mean, lower = lower, upper = upper,  # Numeric columns for the plot
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = font, cex = 1),
      ticks = gpar(cex = 0.88),
      xlab = gpar(cex = 0.88)
    ),
    graph.pos = 2, 
    hrzl_lines = list("2" = gpar(lty = 2),
                      "3" = gpar(lty = 2),
                      "4" = gpar(lty = 2),
                      "5" = gpar(lty = 2),
                      "6" = gpar(lty = 2),
                      "7" = gpar(lty = 2)),
    xlog = TRUE,
    # xticks = log(c(0.80, 1, 1.25)),
    lty.ci = 1,
    col = fpColors(
      box = "maroon4",
      line = "maroon1",
      summary = "magenta4",
      hrz_lines = "gray63"
    ),
    vertices = TRUE,
    xlab = NULL,
    zero = 1
  )

grid.text(
  "Higher outcome risk in no-treatment group < > Higher outcome risk in treatment group",
  x = unit(0.493, "npc"),
  y = unit(0.015, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()

## Subgroup analyses ----
df_results_sg <- df_results_cox %>% 
  filter(grepl("sg_", name)) %>% 
  filter(grepl("mdl_max_adj", model))

df_results_sg$name <- recode(df_results_sg$name,
                              "sg_60orabove" = "Age: 60 years or above",
                              "sg_below60" = "Age: below 60 years",
                              "sg_female" = "Sex: Female",
                              "sg_male" = "Sex: Male",
                              "sg_nonwhite" = "Ethnicity: Non-White",
                              "sg_white" = "Ethnicity: White",
                              "sg_imd1" = "Deprivation: Most deprived",
                              "sg_nonimd1" = "Deprivation: Other",
                              "sg_obese" = "Obese",
                              "sg_nonobese" = "Not Obese",
                              "sg_belowHbA1c59" = "HbA1c: below 59 mmol/mol",
                              "sg_HbA1c59orabove" = "HbA1c: 59 mmol/mol or above"
)
desired_order <- c(
  "Age: 60 years or above",
  "Age: below 60 years",
  "Sex: Female",
  "Sex: Male",
  "Ethnicity: Non-White",
  "Ethnicity: White",
  "Deprivation: Most deprived",
  "Deprivation: Other",
  "Obese",
  "Not Obese",
  "HbA1c: below 59 mmol/mol",
  "HbA1c: 59 mmol/mol or above"
)
df_results_sg <- df_results_sg %>%
  mutate(
    name = factor(name, levels = desired_order)
  ) %>%
  arrange(name)

# build forestplot
base_data <- tibble(
  mean = df_results_sg$hr,
  lower = df_results_sg$conf_low,
  upper = df_results_sg$conf_high,
  name = as.character(df_results_sg$name),
  # events_i = as.character(df_results_sg$sarilumab),
  # events_c = as.character(df_results_sg$tocilizumab),
  estimates = paste0(formatC(df_results_sg$hr, format = "f", digits = 2), 
                     " (", formatC(df_results_sg$conf_low, format = "f", digits = 2), 
                     " - ", formatC(df_results_sg$conf_high, format = "f", digits = 2), ")"))
header <- tibble(
  name = "Subgroup on primary outcome",
  # events_i = "Events Sarilumab\n n (%)",
  # events_c = "Events Tocilizumab\n n (%)",
  estimates = "aHR (95% CI)",
  mean = NA, lower = NA, upper = NA)

fp <- bind_rows(header, base_data)
font <- "sans"

# Create the forest plot
pdf(file = "Fp_results_sg.pdf",
    width = 14,
    height = 12)

fp %>%
  forestplot(
    labeltext = c(name, 
                  # events_i, events_c, 
                  estimates),
    mean = mean, lower = lower, upper = upper,  # Numeric columns for the plot
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = font, cex = 1),
      ticks = gpar(cex = 0.88),
      xlab = gpar(cex = 0.88)
    ),
    graph.pos = 2, 
    hrzl_lines = list("2" = gpar(lty = 2),
                      "4" = gpar(lty = 2),
                      "6" = gpar(lty = 2),
                      "8" = gpar(lty = 2),
                      "10" = gpar(lty = 2),
                      "12" = gpar(lty = 2)),
    xlog = TRUE,
    # xticks = log(c(0.80, 1, 1.25)),
    lty.ci = 1,
    col = fpColors(
      box = "maroon4",
      line = "maroon1",
      summary = "magenta4",
      hrz_lines = "gray63"
    ),
    vertices = TRUE,
    xlab = NULL,
    zero = 1
  )

grid.text(
  "Higher outcome risk in no-treatment group < > Higher outcome risk in treatment group",
  x = unit(0.375, "npc"),
  y = unit(0.010, "npc"),
  gp = gpar(cex = 0.9)
)

dev.off()
