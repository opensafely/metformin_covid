# # # # # # # # # # # # # # # # # # # # #
# This script creates directed acyclic graphs (DAGs) used for the study
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
# install.packages("ggdag")
library(dagitty)
library(ggdag)
library(ggplot2)
library(tidyverse)


# Create the DAG ----
dag <- dagify(covidoutcome ~ metformin + agesexsmoking + ethnicity + rural_urban + imd + cvdcomorb + othercomorb + DMseverity + hsb,
              metformin ~ agesexsmoking + ethnicity + rural_urban + imd + cvdcomorb + othercomorb + DMseverity + hsb,
              DMseverity ~ agesexsmoking,
              cvdcomorb ~ ethnicity + rural_urban + imd + agesexsmoking + hsb,
              othercomorb ~ ethnicity + rural_urban + imd + agesexsmoking + hsb,
              rural_urban ~ imd,
              ethnicity ~ imd,
              hsb ~ ethnicity + rural_urban + imd + agesexsmoking,
              labels = c("covidoutcome" = "COVID outcomes",
                         "metformin" = "Metformin",
                         "agesexsmoking" = "Age, Sex, Smoking",
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
