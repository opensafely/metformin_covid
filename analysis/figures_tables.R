################################################################################
## This script does the following:
# 1. Import released midpoint rounded output data
# 2. Create different descriptive figures 
# 3. Save all output
################################################################################

################################################################################
# Import libraries + functions
################################################################################
library('tidyverse')
library('here')
library('DiagrammeR')
library('ggplot2')

################################################################################
# Create directories for output
################################################################################
fs::dir_create(here::here("output", "data_description"))

################################################################################
# Import the data
################################################################################
df_eligiblity <- read_csv(here("output", "data_description", "n_elig_excluded_midpoint6.csv"))
df_qa <- read_csv(here("output", "data_description", "n_qa_excluded_midpoint6.csv"))
df_completeness <- read_csv(here("output", "data_description", "n_completeness_excluded_midpoint6.csv"))
df_patterns <- read_csv(here("output", "data_description", "n_exp_out_midpoint6.csv"))
df_tbl1 <- read_csv(here("output", "data_description", "table1_midpoint6.csv"))

################################################################################
# Quality assurance flowchart
################################################################################
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

################################################################################
# Completeness criteria flowchart
################################################################################
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

################################################################################
# Eligibility flowchart
################################################################################
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

################################################################################
# Treatment/exposure patterns
################################################################################
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
# Plot 1
plot_1 <- ggplot(df_patterns1, aes(x = reorder(Variable, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Value), hjust = -0.2, size = 5) + 
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal()

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
plot_2 <- ggplot(df_patterns2, aes(x = reorder(Variable, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Value), hjust = -0.2, size = 5) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal()

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
plot_3 <- ggplot(df_patterns3, aes(x = reorder(Variable, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Value), hjust = -0.2, size = 5) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal()

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
                         "Among those without metformin mono 6m, any antidiabetic (+/- metformin) or nothing",
                         "No metformin mono or any other antidiabetic (+/- metformin) within 6m"))
# Plot 4
plot_4 <- ggplot(df_patterns4, aes(x = reorder(Variable, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Value), hjust = -0.2, size = 5) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal()

# Selection for plot 5
df_patterns5 <- df_patterns %>%
  filter(Variable %in% c("COVID hosp or death (after pandemic start)",
                         "COVID hosp or death (after baseline)",
                         "COVID hosp (after baseline)",
                         "COVID death (after baseline)",
                         "COVID diagnosis, pos test or hosp (after baseline)"))
# Plot 5
plot_5 <- ggplot(df_patterns5, aes(x = reorder(Variable, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Value), hjust = -0.2, size = 5) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal()

# Selection for plot 6
df_patterns6 <- df_patterns %>%
  filter(Variable %in% c("Deaths between landmark and pandemic start",
                         "LTFU between landmark and pandemic start",
                         "Metformin (combo) in 6m prior to pandemic start",
                         "Metformin mono in 6m prior to pandemic start"))
# Plot 6
plot_6 <- ggplot(df_patterns6, aes(x = reorder(Variable, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Value), hjust = -0.2, size = 5) +  
  coord_flip() +
  labs(title = "Among all eligible (T2DM)",
       x = "",
       y = "Proportion") +
  theme_minimal()

################################################################################
# Baseline table, based on direct output from table1.R
################################################################################
# df_tbl1

################################################################################
# Save output
################################################################################
# ggsave(filename = here::here("output", "data_description", "plot_6.png"), plot_6, width = 20, height = 20, units = "cm")

