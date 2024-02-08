# NOTES -------------------------------------------------------------------
# Major re-write on 2/25/2023 for new system of trials


# Load Packages -----------------------------------------------------------
# data loading
library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(glue); library(lubridate); library(broom)

# analysis & visualization
library(psycho); library(ggplot2); library(nortest); library(hrbrthemes); library(gtools); 
library(FSA); library(nlstools); library(nlraa); library(directlabels)
# FSA provides an SE calc


# # Clear workspace ---------------------------------------------------------
# rm(list = ls())


# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# Location of the R scripts
code_folder = "Y:/GitHub/Data Analysis/TTS"

# Explict location to save files to:
save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/TTS Graphs"

# Sensitivity cutoff for determining hearing thresholds
TH_cutoff = 1.5

# FA cutoff
FA_cutoff = .41

# Working directory -------------------------------------------------------
setwd(code_folder)


# Groupings ---------------------------------------------------------------

# TTS Pilot Group - 50ms data only b/c mostly mixed
# 2 rats "Green2" (16) & "Orange1" (21) have incomplete data set
Group_TTS_pilot = c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)

# Temporal Integration Pilot Group - all 3 time points and ABRs pre & post HL
# 1 rat "Green12" (191) had threshold shift
Group_1 = c(190, 191, 192, 193)

# Temporal Integration Pilot Group - all 3 time points and ABRs pre & post HL
# 1 rat "Blue2" (266) died during baseline
# 1 rat "Blue3" (267) never got to FA criterion and was dropped
Group_2 = c(265, 266, 267, 268)

# Temporal Integration Pilot Group - all 3 time points and ABRs pre & post HL
# in progress
Group_3 = c(338, 339, 340, 341, 342, 343)



# Import Data -----------------------------------------
# errors, post ABRs and 'maintenance' days are automatically removed
cat("Loading data...")
source(glue("{code_folder}/TTS_data.R"))
cat("done\n")


# # Analysis ----------------------------------------------------------------
# # Graph Hit, FA and Trial Count from summary data
# # calculate hearing threshold (from all data) and remove any trials below hearing level.
# # ISSUE: Plots don't show.
#
# source(glue("{code_folder}/TTS_analysis.R"))
#
# Graphing ----------------------------------------------------------------

cat("Graphing...")

# source(glue("{code_folder}/TTS_graphs.R"))
source(glue("{code_folder}/Temporal Integration.R"))

cat("done\n")

