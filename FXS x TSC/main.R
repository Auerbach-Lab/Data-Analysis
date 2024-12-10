
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



# Load Necessary Datasets -------------------------------------------------
cat("Loading data.\n")
load(paste0(projects_folder, "run_archive.Rdata"), .GlobalEnv)
rat_decoder = fread(glue("{projects_folder}/rat_archive.csv"),
                    select = c("Rat_ID", "DOB", "Sex", "Genotype", "HL_date"))

# Remove bad data
#TODO: deal with multiple runs in a day
dataset = run_archive %>% 
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a dataframe 
  filter(str_detect(warnings_list, pattern = "wrong delay window", negate = TRUE))

# Get core data -----------------------------------------------------------
core_columns = c("date", "rat_name", "rat_ID", "genotype", "sex",
                 "file_name", "experiment", "phase", "task", "detail", 
                 "stim_type", "analysis_type", "block_size", "complete_block_count", 
                 "dprime", "reaction", "FA_percent")

core_data = dataset %>% 
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>% 
  # Only keep relevant Experiments
  filter(experiment %in% c("FXS x TSC")) %>%
  left_join(rat_decoder, by = join_by(rat_ID == Rat_ID)) %>%
  rename(sex = Sex) %>%
  mutate(FXS = str_extract(Genotype, pattern = "Fmr1_[:alpha:]+") %>% str_remove("Fmr1_"),
         TSC = str_extract(Genotype, pattern = "Tsc2_[:alpha:]+") %>% str_remove("Tsc2_"),
         genotype = case_when(FXS == "WT" & TSC == "WT" ~ "Wild-type",
                              FXS == "KO" & TSC == "WT" ~ "FXS only",
                              FXS == "WT" & TSC == "Het" ~ "TSC only",
                              FXS == "KO" & TSC == "Het" ~ "Double KO")) %>%
  unnest_wider(stats) %>%
  select(all_of(core_columns))
