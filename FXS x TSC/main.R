
# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# where to save graphs
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/FXS x TSC"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 1.5

# FA cutoff
FA_cutoff = .41


# Load Packages -----------------------------------------------------------

# data loading
library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(broom); library(glue); library(lubridate)

# analysis & visualization
library(psycho); library(FSA); library(nortest); library(drda); 
library(ggplot2); library(hrbrthemes); library(gtools); library(directlabels)