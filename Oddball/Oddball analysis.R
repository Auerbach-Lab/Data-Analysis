
# Variables ---------------------------------------------------------------
# Location of the datasets
projects_folder = "Z:/Behavior-autoanalysis/"

# FA cutoff
FA_cutoff = .41


# Load Packages -----------------------------------------------------------

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(glue); library(data.table)

# analysis & visualization
library(ggplot2); library(nortest); library(FSA); library(ggbeeswarm)



# Load Necessary Datasets -------------------------------------------------
cat("Loading data...")
load(paste0(projects_folder, "run_archive.Rdata"), .GlobalEnv)
rat_decoder = fread(glue("{projects_folder}/rat_archive.csv"),
                    select = c("Rat_ID", "DOB", "Sex", "Genotype", "HL_date"))
cat("done.\n")

# Individual Trial Data
trials_data = fread(paste0(projects_folder, "Oddball_archive.csv.gz"))

cat("Filtering...")

# Remove bad data
#TODO: deal with multiple runs in a day
dataset = run_archive %>% 
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a dataframe 
  filter(str_detect(warnings_list, pattern = "wrong delay window", negate = TRUE))


# Get core data -----------------------------------------------------------

core_data = dataset %>% 
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>% unnest_wider(stats) %>%
  # Only keep relevant Experiments
  filter(experiment %in% c("Oddball")) %>%
  mutate(task = str_replace_all(task, "CNO \\(3mg/kg\\)", "CNO 3mg/kg")) %>%
  mutate(go = str_extract(file_name, pattern = "^[:digit:]+") %>% as.numeric(), # Get go frequency by extracting 1st (^) number ([:digit:]+) from file_name
         no_go = str_extract(file_name, pattern = "(?<=_)(BBN|[:digit:]+kHz)(?=_)")) %>% # Get No go frequency by extracting either BBN or the 2nd kHz ([:digit:]+) from file_name
  mutate(Freq_pair = case_when(phase == "Tone-BBN" ~ paste(task, go),
                               phase == "Tone-Tone" ~ glue("{task} {go}-{no_go}"),
                               .default = glue("unknown phase: {task} {go}-{no_go}"))) %>%
  left_join(rat_decoder, by = c("rat_ID" = "Rat_ID"))

Tone_on_Tone_rats = filter(core_data, phase == "Tone-Tone")$rat_ID %>% unique

Trials_table = 
  trials_data %>%
  as_tibble() %>%
  filter(UUID %in% core_data$UUID) %>%
  left_join(
    select(core_data, date, rat_name, rat_ID, Genotype,
           UUID, invalid, analysis_type, 
           file_name, experiment, phase, task, detail,
           omit_list, block_size), by = join_by(UUID)) %>%
  # filter out omitted trials
  group_by(UUID) %>% 
  do(filter(., ! Trial_number %in% .$omit_list)) %>%
  ungroup

# Descriptive Stats -------------------------------------------------------

cat("summerizing...")

Summary_table = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  # # Omit days with > 45% FA, i.e. guessing
  # filter(FA_percent < FA_cutoff) %>%
  reframe(trials = mean(trial_count, na.rm = T), 
          blocks = mean(complete_block_count, na.rm = T),
          hit = mean(hit_percent, na.rm = T) * 100, 
          FA = mean(FA_percent, na.rm = T) * 100,
          .by = c(rat_ID, rat_name, Genotype, task, phase, go, no_go))

## Trials calculations ----

Trial_run_summary_table = 
  Trials_table %>%
  as_tibble() %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  reframe(trials = n(), 
          Rxn = mean(`Reaction_(s)`, na.rm = TRUE) * 1000,
          delay = mean(`Delay (s)`, na.rm = TRUE),
          .by = c(UUID, date, rat_ID, rat_name, Genotype, 
                  task, phase, detail, 
                  Stim_ID, Trial_type, `Inten (dB)`, Response)) %>% 
  arrange(UUID, Stim_ID)

multiple_runs_in_a_day = 
  reframe(core_data,
          UUID_count = length(unique(UUID)),
          UUIDs = unique(UUID),
          .by = c(date, rat_ID, rat_name, Genotype, task, phase, detail)) %>%
  filter(UUID_count > 1)

Daily_summary_from_trials = 
  Trials_table %>%
  as_tibble() %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  reframe(trials = n(), 
          Rxn = mean(`Reaction_(s)`, na.rm = TRUE) * 1000,
          delay = mean(`Delay (s)`, na.rm = TRUE),
          .by = c(date, rat_ID, rat_name, Genotype, 
                  task, phase, detail, 
                  Stim_ID, Trial_type, `Inten (dB)`, Response)) %>%
  mutate(total_trials = sum(trials),
         .by = c(date, rat_ID, rat_name, Genotype, task, phase, detail))

# Rat_trial_summary = 
  
  

# Tone-Tone Summary -------------------------------------------------------
Summary_table %>%
  filter(rat_ID > 300) %>%
  filter(task == "Base case") %>%
  reframe(trials = mean(trials, na.rm = T), 
          blocks = mean(blocks, na.rm = T),
          hit = mean(hit, na.rm = T) %>% round(digits = 0), 
          FA = mean(FA, na.rm = T) %>% round(digits = 0),
          .by = c(task, phase, go, no_go)) %>%
  View

# Probe FA count ----------------------------------------------------------

probe_days = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(task %in% c("Probe trials"))

Probe_table = 
  trials_data %>%
  as_tibble() %>%
  # get only days with probes
  filter(UUID %in% probe_days$UUID) %>%
  # get only probe trials (Trial_type = 0 OR Type = 0 OR `Inten (dB)` = 0)
  filter(Trial_type == 0) %>%
  left_join(probe_days, by = join_by(UUID))

Probe_results = 
  Probe_table %>%
  count(UUID, Response) %>% # FAs and CRs are going to be in different rows
  spread(Response, n) %>%   # move them to columns
  left_join(probe_days, by = join_by(UUID)) %>%
  # get table
  reframe(FA_count = sum(FA, na.rm = T),
          Probe_count = sum(CR, FA, na.rm = T),
          Freq = unique(go) %>% sort() %>% str_flatten_comma(),
          .by = c(rat_ID, rat_name, Genotype, task, phase, detail)) %>%
  arrange(rat_ID)

filter(Probe_results, rat_ID > 300) %>% View

## FAs on non-probe trials ---
Probe_FA_count = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(task %in% c("Probe trials")) %>%
  unnest(FA_detailed) %>%
  # # Omit days with > 45% FA, i.e. guessing
  # filter(FA_percent < FA_cutoff) %>%
  reframe(FA_count = sum(FA, na.rm = T),
          Probe_count = sum(trials, na.rm = T),
          Freq = unique(go) %>% sort() %>% str_flatten_comma(),
          .by = c(rat_ID, rat_name, Genotype, task, phase, detail, position)) %>%
  arrange(rat_ID, position)

Probe_FA_count_by_go = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(task %in% c("Probe trials")) %>%
  unnest(FA_detailed) %>%
  # # Omit days with > 45% FA, i.e. guessing
  # filter(FA_percent < FA_cutoff) %>%
  reframe(FA_count = sum(FA, na.rm = T),
          Probe_count = sum(trials, na.rm = T),
          # Freq = unique(go) %>% sort() %>% str_flatten_comma(),
          .by = c(rat_ID, rat_name, Genotype, task, phase, detail, go, position)) %>%
  arrange(rat_ID, go, position)

# filter(Probe_FA_count, rat_ID > 300) %>% View


# Get Reaction times by rat -----------------------------------------------

cat("getting reaction times...")


core_columns = c("date", "rat_name", "rat_ID", "Sex", "Genotype", 
                 "file_name", "experiment", "phase", "task", "detail", 
                 "stim_type", "analysis_type", 
                 "reaction", "FA_percent", "go", "no_go")

Rxn = core_data %>%
  select(all_of(core_columns)) %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  # filter(FA_percent < FA_cutoff) %>%
  # Get Reaction times:
  unnest(reaction) %>%
  rename(Position = `Inten (dB)`)

Rxn_table = Rxn %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail, phase, Position) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>% 
  unique() %>%
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail, phase) %>%
  do(mutate(., Rxn_norm = Rxn/filter(., Position == min(Position))$Rxn,
            Rxn_diff = Rxn - filter(., Position == min(Position))$Rxn)) %>%
  ungroup


Rxn_Tone_On_Tone = Rxn %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail, phase, go, no_go, Position) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>% 
  unique()  %>%
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail, phase, go, no_go) %>%
  do(mutate(., Rxn_norm = Rxn/filter(., Position == min(Position))$Rxn,
            Rxn_diff = Rxn - filter(., Position == min(Position))$Rxn)) %>%
  ungroup


# Rxn analysis ------------------------------------------------------------

# Normality testing
ad.test(Rxn_table$Rxn)

Model_data = Rxn_table %>% 
  filter(task %in% c("Rotating", "CNO 3mg/kg")) %>%
  filter(rat_ID != 100)

## Overall Model
Rxn_overall_model = aov(Rxn ~ task * Position, 
                data = Model_data)

summary(Rxn_overall_model)

kruskal.test(Rxn ~ task, 
             data = Model_data)

Rxn_overall_model_postHoc <- 
  FSA::dunnTest(Rxn ~ interaction(task), 
                data = Rxn_table,
                method = "bonf")

Rxn_overall_model_postHoc$res %>% 
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         P.adj = round(P.adj, digits = 3))



