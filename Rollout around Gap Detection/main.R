# NOTES -------------------------------------------------------------------
# Major re-write on 2/25/2023 for new system of trials


# Load Packages -----------------------------------------------------------
# data loading
library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(glue); library(lubridate); library(broom)

# analysis & visualization
library(psycho); library(ggplot2); library(hrbrthemes);
# FSA provides an SE calc
library(FSA);


# # Clear workspace ---------------------------------------------------------
# rm(list = ls())


# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# Location of the R scripts
code_folder = "Y:/GitHub/Data Analysis/Rollout around Gap Detection"

# Explict location to save files to:
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Rollout Issues"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 1.5

# FA cutoff
FA_cutoff = .41

# Working directory -------------------------------------------------------
setwd(code_folder)


# Import Data -----------------------------------------
# Load Necessary Datasets -------------------------------------------------
load(glue("{projects_folder}/run_archive.Rdata"), .GlobalEnv)


# Analysis ----------------------------------------------------------------
source(glue("{code_folder}/behavior system oddity.R"))
