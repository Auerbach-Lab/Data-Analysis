# NOTES -------------------------------------------------------------------


# Load Packages -----------------------------------------------------------
# data loading
library(data.table); library(openxlsx)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(glue); library(lubridate); library(broom)

# analysis & visualization
library(ggplot2); library(psycho); library(drda); library(hrbrthemes); library(gtools); 
library(FSA); library(directlabels)
# FSA provides an SE calc

# data saving
library(svglite)

# # Clear workspace ---------------------------------------------------------
# rm(list = ls())


# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# Location of excel spreadsheets
excel_folder = "Z:/Project Excels/"

# Location of the R scripts
code_folder = "Y:/GitHub/Data Analysis/Fmr1 SD Octave/"

# Explict location to save files to:
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Fmr1 SD Octave/"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 2.0

# ABR directory
ABR_data_folder = "C:/Users/Noelle/Box/ABR recordings/ABR Results/Fmr1 SD Rats/"

# Box Folder
Box_folder = "C:/Users/Noelle/Box"

# Working directory -------------------------------------------------------
setwd(code_folder)


# Import Data -----------------------------------------
# errors, post ABRs and 'maintenance' days are automatically removed
cat("Loading data...")
source(glue("{code_folder}/data.R"))
cat("done\n")

