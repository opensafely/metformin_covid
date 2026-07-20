# # # # # # # # # # # # # # # # # # # # #
# This script creates directed acyclic graphs (DAGs) used for the study
# Conceptual figure only - does not touch any patient-level data.
# Not a project.yaml action: `ggdag` is not available in the OpenSAFELY r:v2
# image, so run this directly with Rscript/RStudio on your own machine.
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library(dagitty)
library(V8)
library(ggdag)
library(ggplot2)
library(tidyverse)
library(here)

options(dagitty.use.v8 = FALSE)

# Create output directory ----
fs::dir_create(here::here("output", "dag"))
# Create the DAG ----
dag <- dagify(covidoutcome ~ metformin + agesexsmokingbmi + ethnicity + rural_urban + imd + cvdcomorb + othercomorb + DMseverity + hsb,
              metformin ~ agesexsmokingbmi + ethnicity + rural_urban + imd + cvdcomorb + othercomorb + DMseverity + hsb,
              DMseverity ~ agesexsmokingbmi,
              cvdcomorb ~ ethnicity + rural_urban + imd + agesexsmokingbmi + hsb,
              othercomorb ~ ethnicity + rural_urban + imd + agesexsmokingbmi + hsb,
              rural_urban ~ imd,
              ethnicity ~ imd,
              hsb ~ ethnicity + rural_urban + imd + agesexsmokingbmi,
              labels = c("covidoutcome" = "COVID outcomes",
                         "metformin" = "Metformin",
                         "agesexsmokingbmi" = "Age, Sex, Smoking, BMI",
                         "rural_urban" = "Rural/urban",
                         "ethnicity" = "Ethnicity",
                         "cvdcomorb" = "CVDs other\n than T2DM*",
                         "othercomorb" = "Non-CVD comorbidities**",
                         "DMseverity" = "Diabetes severity****",
                         "imd" = "Deprivation",
                         "hsb" = "Health seeking behaviour***"
              ),
              exposure = "metformin",
              outcome = "covidoutcome"
)
plot_dag <- ggdag(dag, text = F, use_labels = "label") +
  theme_void() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )
plot_adjustment_set <- ggdag_adjustment_set(dag, text = F, use_labels = "label") +
  theme_void() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# Save plots ----
ggsave(here::here("output", "dag", "dag_plot.png"), plot_dag, width = 20, height = 16, units = "cm")
ggsave(here::here("output", "dag", "dag_adjustment_set_plot.png"), plot_adjustment_set, width = 20, height = 16, units = "cm")
