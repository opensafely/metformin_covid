# # # # # # # # # # # # # # # # # # # # #
# This script creates directed acyclic graphs (DAGs) used for the study
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library(dagitty)
library(V8)
library(ggdag)
library(ggplot2)
library(tidyverse)

options(dagitty.use.v8 = FALSE)
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
ggdag(dag, text = F, use_labels = "label") +
  theme_void() +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank()
  )
ggdag_adjustment_set(dag, text = F, use_labels = "label") +
  theme_void() +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank()
  )
