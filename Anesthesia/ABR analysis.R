
# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# where to save graphs
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Anesthesia"

# Load Packages -----------------------------------------------------------

# data loading
library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(broom); library(glue); library(lubridate)

# analysis & visualization
library(psycho); library(FSA); library(nortest); library(drda); 
library(ggplot2); library(hrbrthemes);


# Load Necessary Datasets -------------------------------------------------
cat("Loading data.\n")
data = fread("RMS_diff_output.csv")


# MANOVA ------------------------------------------------------------------

aov = aov(RMS_diff ~ Intensity * Drug_Assigned, data = data)

summary(aov)

broom::tidy(TukeyHSD(aov, which = "Drug_Assigned")) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  filter(sig != " ")


# post_hoc = 
  data %>%
  group_by(Intensity) %>%
    do(
      TukeyHSD(aov(RMS_diff ~ Drug_Assigned, data = .)) %>% tidy
    ) %>%
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
    filter(sig != " ")
  
  