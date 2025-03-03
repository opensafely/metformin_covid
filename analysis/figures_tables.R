####
## This script does the following:
# 1. Import released midpoint rounded disclosure-safe output data
# 2. Create different descriptive figures 
# 3. Save all output
####

# Import libraries ----
library('dplyr')
library('tidyr')
library('readr')
library('stringr')
library('here')
library('DiagrammeR')
library('ggplot2')
library('RColorBrewer')
library('gt')

####
# Create directories for output ----
####
fs::dir_create(here::here("output", "data_description"))

####
# Import the data ----
####
df_eligiblity <- read_csv(here("output", "data_release_20250303", "n_elig_excluded_midpoint6.csv"))
df_qa <- read_csv(here("output", "data_release_20250303", "n_qa_excluded_midpoint6.csv"))
df_completeness <- read_csv(here("output", "data_release_20250303", "n_completeness_excluded_midpoint6.csv"))
df_patterns <- read_csv(here("output", "data_release_20250303", "n_exp_out_midpoint6.csv"))
df_tbl1 <- read_csv(here("output", "data_release_20250303", "table1_midpoint6.csv"))

####
# Quality assurance flowchart ----
####
# Create node labels
node_labels <- paste(df_qa$Variable, '\n', df_qa$Value)

# Create flowchart
graph_qa <- grViz(
  sprintf(
    "digraph flowchart {
      graph [layout = dot, rankdir = TB]
      node [shape = box, style = filled, fontname = Arial, fontsize = 14]

      A [label = '%s', fillcolor = lightskyblue]
      B [label = 'Exclusions\n%s\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
      C [label = '%s', fillcolor = lightskyblue]

      A -> B
      B -> C
    }",
    node_labels[1], node_labels[2],
    node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
    node_labels[8], node_labels[9]
  )
)

####
# Completeness criteria flowchart ----
####
# Create node labels
node_labels <- paste(df_completeness$Variable, '\n', df_completeness$Value)

# Create flowchart
graph_completeness <- grViz(
  sprintf(
    "digraph flowchart {
      graph [layout = dot, rankdir = TB]
      node [shape = box, style = filled, fontname = Arial, fontsize = 14]

      A [label = '%s', fillcolor = lightskyblue]
      B [label = 'Exclusions\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
      C [label = '%s', fillcolor = lightskyblue]

      A -> B
      B -> C
    }",
    node_labels[1], node_labels[2],
    node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
    node_labels[8]
  )
)

####
# Eligibility flowchart ----
####
# Create node labels
node_labels <- paste(df_eligiblity$Variable, '\n', df_eligiblity$Value)

# Create flowchart
graph_elig <- grViz(
  sprintf(
    "digraph flowchart {
      graph [layout = dot, rankdir = TB]
      node [shape = box, style = filled, fontname = Arial, fontsize = 14]

      A [label = '%s', fillcolor = lightskyblue]
      B [label = '%s', fillcolor = lightyellow]
      C [label = 'Exclusions 1\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
      D [label = 'Exclusions 2: Events between baseline and landmark\n%s\n%s\n%s\n%s\n%s\n%s', fillcolor = lightyellow]
      E [label = '%s', fillcolor = lightskyblue]

      A -> B
      B -> C
      C -> D
      D -> E
    }",
    node_labels[1], node_labels[2], 
    node_labels[3], node_labels[4], node_labels[5], node_labels[6], node_labels[7],
    node_labels[8], node_labels[9], node_labels[10], node_labels[11], node_labels[12], 
    node_labels[13], node_labels[14], node_labels[15], 
    node_labels[16], node_labels[17], node_labels[18], node_labels[19], node_labels[20], node_labels[21],
    node_labels[22]
  )
)

####
# Treatment/exposure patterns ----
####
# Add a column for denominator = N
N <- df_eligiblity %>%
  filter(Variable == "After applying eligibility criteria") %>%
  select(Value) %>%
  pull()
df_patterns$Proportion <- df_patterns$Value / N

# Selection for plot 1
df_patterns1 <- df_patterns %>%
  filter(Variable %in% c("Metformin (combo) within 6m", "Metformin mono within 6m",
                         "Metformin (combo) within 3m", "Metformin mono within 3m",
                         "Metformin (combo) anytime", "Metformin mono anytime"))
df_patterns1 <- df_patterns1 %>%
  select(-...1) %>% 
  distinct() # take out the duplicate variable, that was released for double-checking

# Plot 1
df_patterns1$Variable <- factor(df_patterns1$Variable, levels = df_patterns1$Variable[order(df_patterns1$Proportion, decreasing = TRUE)])
plot_1 <- ggplot(df_patterns1, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Value), hjust = -0.2, size = 4) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +  # Remove unnecessary legend
  scale_fill_brewer(palette = "Blues")

# Selection for plot 1.1
df_patterns1.1 <- df_patterns %>%
  select(Variable, Value) %>% 
  filter(Variable %in% c("Median time from T2DM diagnosis to metformin (combo) start",
                         "IQR lower bound: T2DM diagnosis to metformin (combo)",
                         "IQR upper bound: T2DM diagnosis to metformin (combo)")) %>% 
  pivot_wider(
    names_from = "Variable",
    values_from = "Value"
  ) %>% 
  rename(Median = "Median time from T2DM diagnosis to metformin (combo) start",
         IQR_lower = "IQR lower bound: T2DM diagnosis to metformin (combo)",
         IQR_upper = "IQR upper bound: T2DM diagnosis to metformin (combo)") %>% 
  mutate(Variable = "Median, IQR: T2DM diagnosis to metformin (combo) start")
# Plot 1.1
plot_1.1 <- ggplot(df_patterns1.1, aes(x = Variable, y = Median)) +
  geom_errorbar(aes(ymin = IQR_lower, ymax = IQR_upper), width = 0.2, color = "red") +
  geom_point(aes(y = Median), color = "black", size = 3) +
  labs(title = "Boxplot: T2DM Diagnosis to Metformin Start",
       x = "",
       y = "Days") +
  theme_minimal()

# Selection for plot 2
df_patterns2 <- df_patterns %>%
  filter(Variable %in% c("Among those without metformin (combo) anytime, DPP4",
                         "Among those without metformin (combo) anytime, TZD",
                         "Among those without metformin (combo) anytime, SGLT2",
                         "Among those without metformin (combo) anytime, sulfonylurea",
                         "Among those without metformin (combo) anytime, GLP1",
                         "Among those without metformin (combo) anytime, meglitinide",
                         "Among those without metformin (combo) anytime, alpha-glucosidase",
                         "Among those without metformin (combo) anytime, insulin",
                         "No metformin (combo) or any other antidiabetic anytime"))
# Plot 2
df_patterns2$Variable <- factor(df_patterns2$Variable, levels = df_patterns2$Variable[order(df_patterns2$Proportion, decreasing = TRUE)])
plot_2 <- ggplot(df_patterns2, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Value), hjust = -0.2, size = 3) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +  # Remove unnecessary legend
  scale_fill_brewer(palette = "Blues")


# Selection for plot 3
df_patterns3 <- df_patterns %>%
  filter(Variable %in% c("Among those without metformin (combo) 6m, DPP4",
                         "Among those without metformin (combo) 6m, TZD",
                         "Among those without metformin (combo) 6m, SGLT2",
                         "Among those without metformin (combo) 6m, sulfonylurea",
                         "Among those without metformin (combo) 6m, GLP1",
                         "Among those without metformin (combo) 6m, meglitinide",
                         "Among those without metformin (combo) 6m, alpha-glucosidase",
                         "Among those without metformin (combo) 6m, insulin",
                         "Among those without metformin (combo) 6m, any antidiabetic (mono)",
                         "No metformin (combo) or any other antidiabetic within 6m"))
# Plot 3
df_patterns3$Variable <- factor(df_patterns3$Variable, levels = df_patterns3$Variable[order(df_patterns3$Proportion, decreasing = TRUE)])
plot_3 <- ggplot(df_patterns3, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Value), hjust = -0.2, size = 3) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +  # Remove unnecessary legend
  scale_fill_brewer(palette = "Blues")


# Selection for plot 4
df_patterns4 <- df_patterns %>%
  filter(Variable %in% c("Among those without metformin mono 6m, DPP4 (+/- metformin)",
                         "Among those without metformin mono 6m, TZD (+/- metformin)",
                         "Among those without metformin mono 6m, SGLT2 (+/- metformin)",
                         "Among those without metformin mono 6m, sulfonylurea (+/- metformin)",
                         "Among those without metformin mono 6m, GLP1 (+/- metformin)",
                         "Among those without metformin mono 6m, meglitinide (+/- metformin)",
                         "Among those without metformin mono 6m, alpha-glucosidase (+/- metformin)",
                         "Among those without metformin mono 6m, insulin (+/- metformin)",
                         "Among those without metformin mono 6m, any antidiabetic (+/- metformin)",
                         # "Among those without metformin mono 6m, any antidiabetic (+/- metformin) or nothing",
                         "No metformin mono or any other antidiabetic (+/- metformin) within 6m"))
# Plot 4
df_patterns4$Variable <- factor(df_patterns4$Variable, levels = df_patterns4$Variable[order(df_patterns4$Proportion, decreasing = TRUE)])
plot_4 <- ggplot(df_patterns4, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Value), hjust = -0.2, size = 3) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +  # Remove unnecessary legend
  scale_fill_brewer(palette = "Blues")


# Selection for plot 5
df_patterns5 <- df_patterns %>%
  filter(Variable %in% c("COVID hosp or death (after pandemic start)",
                         "COVID hosp or death (after baseline)",
                         "COVID hosp (after baseline)",
                         "COVID death (after baseline)",
                         "COVID diagnosis, pos test or hosp (after baseline)"))
# Plot 5
df_patterns5$Variable <- factor(df_patterns5$Variable, levels = df_patterns5$Variable[order(df_patterns5$Proportion, decreasing = TRUE)])
plot_5 <- ggplot(df_patterns5, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Value), hjust = -0.2, size = 3) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +  # Remove unnecessary legend
  scale_fill_brewer(palette = "Blues")


# Selection for plot 6
df_patterns6 <- df_patterns %>%
  filter(Variable %in% c("Deaths between landmark and pandemic start",
                         "LTFU between landmark and pandemic start",
                         "Metformin (combo) in 6m prior to pandemic start",
                         "Metformin mono in 6m prior to pandemic start"))
# Plot 6
df_patterns6$Variable <- factor(df_patterns6$Variable, levels = df_patterns6$Variable[order(df_patterns6$Proportion, decreasing = TRUE)])
plot_6 <- ggplot(df_patterns6, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = Value), hjust = -0.2, size = 3) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +  # Remove unnecessary legend
  scale_fill_brewer(palette = "Blues")


####
# Baseline table, based on direct output from table1.R ----
####
# Step 1: Format the stats and create one label
df_table1 <- df_tbl1 %>%
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

# Step 2: Create 'Unknown' level for each variable and each level
unknown_levels <- df_table1 %>%
  group_by(var_label, by) %>%
  summarise(variable_levels = "Unknown", stat_label = as.character(first(N_miss)), .groups = "drop")
df_table1 <- df_table1 %>%
  bind_rows(unknown_levels) %>%
  arrange(var_label, variable_levels, by)

# Step 3: Reshape the data to wide format
baseline_table <- df_table1 %>% 
  select(var_label, variable_levels, by, stat_label) %>% 
  distinct() %>%  # Ensure unique rows before pivoting
  pivot_wider(
    names_from = by, 
    values_from = stat_label,
    values_fn = list(stat_label = first)  # Take the first value if duplicates exist
  )

# Step 4: Arrange and clean the output
custom_order <- c("Total N", "Age", "Age groups", "Sex", "Ethnicity", "Deprivation", "Region", "Rural/urban", "Smoking status", "Care/nursing home resident"
                  # , "Body Mass Index > 40 kg/m^2", "HbA1c in mmol/mol", "TC/Chol ratio",
                  , "HbA1c categories in mmol/mol", "TC/HDL ratio categories", 
                  "History of acute myocardial infarct", "History of stroke", "History of other arterial embolism", "History of venous thromboembolism",
                  "History of heart failure", "History of angina pectoris", "History of dementia", "History of cancer", "History of arterial hypertension",
                  "History of depression", "History of COPD", "History of liver disease", "History of CKD", "History of PCOS", "History of prediabetes",
                  "Diabetes complication", "Healthcare worker",
                  "Any metformin prescription within 6m prior to pandemic start", 
                  "Starting metformin anytime (after landmark, respectively; incl. combo)",
                  "COVID hosp or death", "COVID hosp", "COVID death", "Any covid diagnosis, pos test or hosp",
                  "Death between landmark and pandemic start",
                  "LTFU between landmark and pandemic start")

# Order the Variables
baseline_table <- baseline_table %>%
  mutate(var_label = factor(var_label, 
                            levels = c(custom_order, setdiff(unique(var_label), custom_order)), 
                            ordered = TRUE)) %>%
  arrange(var_label) # Unspecified variables are kept but pushed to the end

baseline_table <- baseline_table %>%
  # Replace NAs in variable_levels with a placeholder for correct ordering re "Unknown"
  mutate(variable_levels = ifelse(is.na(variable_levels), "NA", variable_levels)) %>%
  # Move 'Unknown' to the bottom of 'variable_levels' within each 'var_label'
  mutate(variable_levels = factor(variable_levels, levels = c(setdiff(unique(variable_levels), "Unknown"), "Unknown"))) %>% 
  # Now sort
  arrange(var_label, variable_levels) %>%
  ungroup()

# Drop "Unknown" rows where var_label has already Unknown as a level or not needed
baseline_table <- baseline_table %>%
  filter(!(var_label == "Deprivation" & variable_levels == "Unknown")) %>% 
  filter(!(var_label == "Total N" & variable_levels == "Unknown")) %>% 
  filter(!(var_label == "Care/nursing home resident" & variable_levels == "Unknown")) %>% 
  filter(!(grepl("History of", var_label) & variable_levels == "Unknown")) %>% 
  filter(!(var_label == "Diabetes complication" & variable_levels == "Unknown")) %>% 
  filter(!(var_label == "Healthcare worker" & variable_levels == "Unknown")) %>% 
  filter(!(grepl("metformin", var_label) & variable_levels == "Unknown")) %>% 
  filter(!(grepl("COVID", var_label) & variable_levels == "Unknown")) %>% 
  filter(!(grepl("Any covid diagnosis", var_label) & variable_levels == "Unknown")) %>% 
  filter(!(grepl("Death between", var_label) & variable_levels == "Unknown")) %>% 
  filter(!(grepl("LTFU between", var_label) & variable_levels == "Unknown"))
  # filter(!(var_label == "Body Mass Index > 40 kg/m^2" & variable_levels == "Unknown"))
  
# Only keep the variable label for the first row of each variable group
baseline_table <- baseline_table %>%
  mutate(
    var_label = if_else(
      row_number() == 1 | var_label != lag(var_label),
      var_label,
      NA_character_
    )
  )

# some further cosmetics: Hide NAs
baseline_table <- baseline_table %>%
  # Hide NAs
  mutate(variable_levels = recode(variable_levels, "NA" = "")) %>% 
  mutate(var_label = replace_na(var_label, ""))


## Create the table
table_gt <- baseline_table %>%
  select(var_label, variable_levels, `Metformin mono`, Nothing, `Other antidiabetic (incl. combo with metformin)`) %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  cols_label(
    var_label = "Characteristic",
    variable_levels = "",  # Empty label
    `Metformin mono` = "Metfin mono",
    Nothing = "Nothing",
    `Other antidiabetic (incl. combo with metformin)` = "Other OAD (incl. metfin combo)"
  ) %>%
  tab_options(
    table.font.size = px(12)
  )


####
# Save output ----
####
# ggsave(filename = here::here("output", "data_description", "plot_6.png"), plot_6, width = 20, height = 20, units = "cm")
