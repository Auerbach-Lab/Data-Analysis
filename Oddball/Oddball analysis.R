
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
library(ggplot2); library(nortest); library(FSA); library(ggbeeswarm); library(effectsize)

# Functions ---------------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}

# Load Necessary Datasets -------------------------------------------------
cat("Loading data...")
load(paste0(projects_folder, "run_archive.Rdata"), .GlobalEnv)
rat_decoder = fread(glue("{projects_folder}/rat_archive.csv"),
                    select = c("Rat_ID", "DOB", "Sex", "Genotype", "HL_date"))

# Individual Trial Data
trials_data = fread(paste0(projects_folder, "Oddball_archive.csv.gz"))
cat("done.\n")

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


## Tone on Tone group ----
Tone_on_Tone_rats = filter(core_data, phase == "Tone-Tone")$rat_ID %>% unique

## CNO on FXS Oddball Group 1 ----
CNO_rats = core_data %>%
  filter(str_detect(Genotype, pattern = "Fmr1-LE.*")) %>%
  # only keep relevant trial types/days
  filter(task %in% c("Base case")) %>%
  filter(detail %in% c("CNO 3mg/kg")) %>%
  # Drop BP2 who was a pilot animal
  filter(rat_name != "BP2") %>%
  # # only keep rats that have been treated with CNO in FXS Oddball Group 1 by name
  # filter(rat_name %in% c("Green1", "Green2", "Green3", "Green4", "Lime5", "Purple5")) %>%
  .$rat_ID %>% unique()

CNO_dates = core_data %>% 
  filter(rat_ID %in% CNO_rats & detail == "CNO 3mg/kg") %>% 
  group_by(rat_ID) %>% 
  do(arrange(., date) %>% head(n = 1)) %>% 
  rename(treatment_date = date) %>%
  select(rat_ID, treatment_date)

## Multiple runs in a day; can be filtered out later ----
multiple_runs_in_a_day = 
  reframe(core_data,
          UUID_count = length(unique(UUID)),
          UUIDs = unique(UUID),
          .by = c(date, rat_ID, rat_name, Genotype, task, phase, detail)) %>%
  filter(UUID_count > 1)


# Get Trial Data ----------------------------------------------------------
Trials_table = 
  trials_data %>%
  as_tibble() %>%
  filter(UUID %in% core_data$UUID) %>%
  left_join(
    select(core_data, date, rat_name, rat_ID, Genotype, Sex,
           UUID, invalid, analysis_type, 
           file_name, experiment, phase, task, detail,
           omit_list, block_size), by = join_by(UUID)) %>%
  # filter out omitted trials
  group_by(UUID) %>% 
  do(filter(., ! Trial_number %in% .$omit_list)) %>%
  ungroup

## Trials calculations; use Daily_summary_from_trials ----
Trial_run_summary_table = 
  Trials_table %>%
  as_tibble() %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  reframe(trials = n(), 
          Rxn = mean(`Reaction_(s)`, na.rm = TRUE) * 1000,
          delay = mean(`Delay (s)`, na.rm = TRUE),
          .by = c(UUID, date, rat_ID, rat_name, Genotype, Sex,
                  task, phase, detail, 
                  Stim_ID, Trial_type, `Inten (dB)`, Response)) %>% 
  arrange(UUID, Stim_ID)

Daily_summary_from_trials = 
  Trials_table %>%
  as_tibble() %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  mutate(go = str_extract(file_name, pattern = "^[:digit:]+") %>% as.numeric(), # Get go frequency by extracting 1st (^) number ([:digit:]+) from file_name
         no_go = str_extract(file_name, pattern = "(?<=_)(BBN|[:digit:]+kHz)(?=_)")) %>% # Get No go frequency by extracting either BBN or the 2nd kHz ([:digit:]+) from file_name
  reframe(trials = n(), 
          Rxn = mean(`Reaction_(s)`, na.rm = TRUE) * 1000,
          delay = mean(`Delay (s)`, na.rm = TRUE),
          .by = c(date, rat_ID, rat_name, Genotype, Sex,
                  task, phase, detail,
                  Trial_type, Response, go, no_go)) %>%
  mutate(total_trials = sum(trials),
         percent = trials/total_trials,
         .by = c(date, rat_ID, rat_name, Genotype, Sex, task, phase, detail, go, no_go))

## Trials by position table -----
Daily_summary_from_trials_by_position = 
  Trials_table %>%
  as_tibble() %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  mutate(go = str_extract(file_name, pattern = "^[:digit:]+") %>% as.numeric(), # Get go frequency by extracting 1st (^) number ([:digit:]+) from file_name
         no_go = str_extract(file_name, pattern = "(?<=_)(BBN|[:digit:]+kHz)(?=_)")) %>% # Get No go frequency by extracting either BBN or the 2nd kHz ([:digit:]+) from file_name
  reframe(trials = n(), 
          Rxn = mean(`Reaction_(s)`, na.rm = TRUE) * 1000,
          delay = mean(`Delay (s)`, na.rm = TRUE),
          .by = c(date, rat_ID, rat_name, Genotype, Sex, 
                  task, phase, detail,
                  Trial_type, `Inten (dB)`, Response, go, no_go)) %>%
  mutate(total_trials = sum(trials),
         .by = c(date, rat_ID, rat_name, Genotype, Sex, task, phase, detail, go, no_go)) %>%
  mutate(position_trials = sum(trials),
         .by = c(date, rat_ID, rat_name, Genotype, Sex, task, phase, detail, go, no_go, `Inten (dB)`)) %>%
  rename(Position = `Inten (dB)`)

Daily_trial_stats = Daily_summary_from_trials_by_position %>%
  mutate(percent = trials/position_trials,
         .by = c(date, rat_ID, rat_name, Genotype, Sex, task, phase, detail, go, no_go, Position, Response)) %>% 
  pivot_wider(names_from = Response,
              names_glue = "{Response}_percent",
              values_from = percent)

## FA table from runs ----
FA_table = Trials_table %>%
  as_tibble() %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  mutate(go = str_extract(file_name, pattern = "^[:digit:]+") %>% as.numeric(), # Get go frequency by extracting 1st (^) number ([:digit:]+) from file_name
         no_go = str_extract(file_name, pattern = "(?<=_)(BBN|[:digit:]+kHz)(?=_)")) %>% # Get No go frequency by extracting either BBN or the 2nd kHz ([:digit:]+) from file_name
  filter(Response == "FA")


# Descriptive Stats -------------------------------------------------------

cat("summerizing...")

Summary_table_run_data = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  # # Omit days with > 45% FA, i.e. guessing
  # filter(FA_percent < FA_cutoff) %>%
  reframe(trials = mean(trial_count, na.rm = T), 
          blocks = mean(complete_block_count, na.rm = T),
          hit = mean(hit_percent, na.rm = T) * 100, 
          FA = mean(FA_percent, na.rm = T) * 100,
          .by = c(rat_ID, rat_name, Genotype, task, phase, go, no_go))

Summary_table = Daily_trial_stats %>%
  # # Omit days with > 45% FA, i.e. guessing
  # filter(FA_percent < FA_cutoff) %>%
  reframe(trials = mean(total_trials, na.rm = T), 
          hit = mean(Hit_percent, na.rm = T) * 100, 
          FA = mean(FA_percent, na.rm = T) * 100,
          miss = mean(Miss_percent, na.rm = T) * 100,
          .by = c(rat_ID, rat_name, Genotype, task, phase, go, no_go))


# Tone-Tone Summary -------------------------------------------------------
Summary_table_run_data %>%
  filter(rat_ID %in% Tone_on_Tone_rats) %>%
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


# Get Reaction times by rat -----------------------------------------------

cat("getting reaction times...")

Rxn_table = Daily_summary_from_trials_by_position %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail, phase, go, Position) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE)) %>% 
  unique() %>%
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail, phase) %>%
  do(mutate(., Rxn_norm = Rxn/filter(., Position == min(Position))$Rxn,
            Rxn_diff = Rxn - filter(., Position == min(Position))$Rxn)) %>%
  ungroup



# Reaction for FXS baseline -----------------------------------------------
FXS_core_baseline_data = 
  Daily_summary_from_trials %>%
  filter(str_detect(Genotype, pattern = "Fmr1")) %>%
  filter(phase == "Tone-BBN") %>%
  filter(Response == "Hit") %>%
  filter(task == "Base case") %>%
  # select only baselines
  filter(str_detect(detail, pattern = "(4-6|Round [:digit:])"))

FXS_core_baseline_data_trial_count = 
  FXS_core_baseline_data %>%
    # Get Averages for each rat
  reframe(go_freq = unique(go),
          go_count = n(),
          .by = c(rat_ID, rat_name, Sex, Genotype, go)) %>%
  reframe(min_count = min(go_count),
          .by = c(rat_ID, rat_name, Sex, Genotype)) %>%
  select(rat_ID, min_count)
  
 
FXS_baseline_hit_reaction =
  FXS_core_baseline_data %>%
    select(-task, -phase, -detail) %>%
  # add trial count to select down to for each rat
  left_join(FXS_core_baseline_data_trial_count,
            by = join_by(rat_ID)) %>%
  nest_by(date, rat_ID, rat_name, go, min_count) %>%
  # get even amount of trials for each frequency
  group_by(rat_ID, rat_name, go) %>% 
  do(arrange(., desc(date)) %>%
       # Select most recent days for each frequency and task
       head(n = unique(.$min_count)) %>%
       unnest(data)) %>%
  ungroup %>%
  # Get Averages for each rat
  reframe(Rxn = mean(Rxn, na.rm = TRUE),
          n = n(),
          .by = c(rat_ID, rat_name, Sex, Genotype))


# DREADD testing (PFC) ----------------------------------------------------------

# Normality testing
ad.test(Rxn_table$Rxn)

Model_data_drug_piloting = Rxn_table %>% 
  filter(task %in% c("Rotating", "CNO 3mg/kg")) %>%
  filter(rat_ID != 100)

## Overall Model
drug_piloting_Rxn_overall_model = 
  aov(Rxn ~ task * Position,
      data = Model_data_drug_piloting)

summary(drug_piloting_Rxn_overall_model)

kruskal.test(Rxn ~ task, 
             data = Model_data_drug_piloting)

drug_piloting_Rxn_overall_model_postHoc <- 
  FSA::dunnTest(Rxn ~ interaction(task), 
                data = Rxn_table,
                method = "bonf")

drug_piloting_Rxn_overall_model_postHoc$res %>% 
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         P.adj = round(P.adj, digits = 3))


# AC inhibition -----------------------------------------------------------

## TODO: filter to 5 before and after treatment
# left_join(CNO_dates,
#           by = join_by(rat_ID))
#   filter(date < treatment_date)

AC_Model_data = Daily_summary_from_trials %>% 
  filter(rat_ID %in% CNO_rats) %>%
  filter(detail %in% c("Round 1", "Round 2", "Round 3", "Between Treatment", "CNO 3mg/kg")) %>%
  mutate(detail = str_replace(detail, pattern = "Round [:digit:]", replacement = "Baseline")) %>%
  nest_by(date, rat_ID, rat_name, go, task, detail) %>%
  arrange(rat_ID, go, date) %>%
  group_by(rat_ID, rat_name, go, detail) %>% 
  do(arrange(., desc(date)) %>% 
       # Select most recent days for each frequency and task
       head(n = 5) %>%
       unnest(data))  %>%
  ungroup %>%
  mutate(detail = factor(detail, 
                         levels = c("Baseline", "Between Treatment", "CNO 3mg/kg")),
         Genotype = str_remove(Genotype, pattern = "Fmr1-LE_"))

AC_Model_data_by_position = Daily_summary_from_trials_by_position %>% 
  filter(rat_ID %in% CNO_rats) %>%
  filter(detail %in% c("Round 1", "Round 2", "Round 3", "Between Treatment", "CNO 3mg/kg")) %>%
  mutate(detail = str_replace(detail, pattern = "Round [:digit:]", replacement = "Baseline")) %>%
  nest_by(date, rat_ID, rat_name, go, task, detail) %>%
  arrange(rat_ID, go, date) %>%
  group_by(rat_ID, rat_name, go, detail) %>% 
  do(arrange(., desc(date)) %>% 
       # Select most recent days for each frequency and task
       head(n = 5) %>%
       unnest(data))  %>%
  ungroup %>%
  mutate(percent = trials/position_trials,
         .by = c(date, rat_ID, rat_name, Genotype, Sex, 
                 task, detail, go, Position, Response)) %>%
  mutate(detail = factor(detail, 
                         levels = c("Baseline", "Between Treatment", "CNO 3mg/kg")),
         Genotype = str_remove(Genotype, pattern = "Fmr1-LE_"))

AC_FA_data = FA_table %>% 
  filter(rat_ID %in% CNO_rats) %>%
  filter(detail %in% c("Round 1", "Round 2", "Round 3", "Between Treatment", "CNO 3mg/kg")) %>%
  mutate(detail = str_replace(detail, pattern = "Round [:digit:]", replacement = "Baseline")) %>%
  nest_by(date, rat_ID, rat_name, go, task, detail) %>%
  arrange(rat_ID, go, date) %>%
  group_by(rat_ID, rat_name, go, detail) %>% 
  do(arrange(., desc(date)) %>% 
       # Select most recent days for each frequency and task
       head(n = 5) %>%
       unnest(data))  %>%
  ungroup %>%
  mutate(detail = factor(detail, 
                         levels = c("Baseline", "Between Treatment", "CNO 3mg/kg")),
         Genotype = str_remove(Genotype, pattern = "Fmr1-LE_"))
  

# AC_Model_data_summary = 
AC_Model_data %>%
  reframe(days = unique(date) %>% length(),
          .by = c(rat_ID, rat_name, Genotype, Sex, task, detail, go)) %>%
  arrange(rat_ID, detail, go)

## Hit/Miss/FA overall ----
### FA ----
FA_AC_aov = 
  aov(percent ~ detail * Genotype,
    data = AC_Model_data %>%
      filter(task == "Base case") %>%
      mutate(detail = factor(detail, 
                             levels = c("Baseline", "Between Treatment", "CNO 3mg/kg")),
             Genotype = str_remove(Genotype, pattern = "Fmr1-LE_")) %>%
      reframe(percent = mean(percent, na.rm = TRUE),
              .by = c(rat_ID, rat_name, Genotype, Sex, 
                      task, detail, go, Response)) %>%
      filter(Response == "FA"))

Parametric_Check(FA_AC_aov)

summary(FA_AC_aov)



#### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "detail" # Main effects
), 
function(x) kruskal.test(reformulate(x, "percent"),
                         data = AC_Model_data %>%
                           filter(task == "Base case") %>%
                           reframe(percent = mean(percent, na.rm = TRUE),
                                   .by = c(rat_ID, rat_name, Genotype, Sex, 
                                           task, detail, go, Response)) %>%
                           filter(Response == "FA"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

#### Post-Hoc Dunn's Test ----
FSA::dunnTest(percent ~ interaction(Genotype, detail),
              data = AC_Model_data %>%
                filter(task == "Base case") %>%
                reframe(percent = mean(percent, na.rm = TRUE),
                        .by = c(rat_ID, rat_name, Genotype, Sex, 
                                task, detail, go, Response)) %>%
                filter(Response == "FA"),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         detail1 = str_split_fixed(Comp1, '\\.', 2)[,2],
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         detail2 = str_split_fixed(Comp2, '\\.', 2)[,2]) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(geno1 == geno2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2)

### Hit ----
Hit_AC_aov = 
  aov(percent ~ detail * Genotype,
      data = AC_Model_data %>%
        filter(task == "Base case") %>%
        reframe(percent = mean(percent, na.rm = TRUE),
                .by = c(rat_ID, rat_name, Genotype, Sex, 
                        task, detail, go, Response)) %>%
        filter(Response == "Hit"))

Parametric_Check(Hit_AC_aov)

summary(Hit_AC_aov)

broom::tidy(TukeyHSD(Hit_AC_aov)) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  filter(sig != " ") 


#### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "detail" # Main effects
), 
function(x) kruskal.test(reformulate(x, "percent"),
                         data = AC_Model_data %>%
                           filter(task == "Base case") %>%
                           reframe(percent = mean(percent, na.rm = TRUE),
                                   .by = c(rat_ID, rat_name, Genotype, Sex, 
                                           task, detail, go, Response)) %>%
                           filter(Response == "Hit"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

#### Post-Hoc Dunn's Test ----
FSA::dunnTest(percent ~ interaction(Genotype, detail),
              data = AC_Model_data %>%
                filter(task == "Base case") %>%
                reframe(percent = mean(percent, na.rm = TRUE),
                        .by = c(rat_ID, rat_name, Genotype, Sex, 
                                task, detail, go, Response)) %>%
                filter(Response == "Hit"),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         detail1 = str_split_fixed(Comp1, '\\.', 2)[,2],
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         detail2 = str_split_fixed(Comp2, '\\.', 2)[,2]) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(geno1 == geno2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2)

### Miss ----
Miss_AC_aov = 
  aov(percent ~ detail * Genotype,
      data = AC_Model_data %>%
        filter(task == "Base case") %>%
        reframe(percent = mean(percent, na.rm = TRUE),
                .by = c(rat_ID, rat_name, Genotype, Sex, 
                        task, detail, go, Response)) %>%
        filter(Response == "Miss"))

Parametric_Check(Miss_AC_aov)

summary(Miss_AC_aov)

#### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "detail" # Main effects
), 
function(x) kruskal.test(reformulate(x, "percent"),
                         data = AC_Model_data %>%
                           filter(task == "Base case") %>%
                           reframe(percent = mean(percent, na.rm = TRUE),
                                   .by = c(rat_ID, rat_name, Genotype, Sex, 
                                           task, detail, go, Response)) %>%
                           filter(Response == "Miss"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

#### Post-Hoc Dunn's Test ----
FSA::dunnTest(percent ~ interaction(Genotype, detail),
              data = AC_Model_data %>%
                filter(task == "Base case") %>%
                reframe(percent = mean(percent, na.rm = TRUE),
                        .by = c(rat_ID, rat_name, Genotype, Sex, 
                                task, detail, go, Response)) %>%
                filter(Response == "Miss"),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         detail1 = str_split_fixed(Comp1, '\\.', 2)[,2],
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         detail2 = str_split_fixed(Comp2, '\\.', 2)[,2]) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(geno1 == geno2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2)

## Hit by position ----
Hit_position_AC_aov = 
  aov(percent ~ detail * Genotype,
      data = AC_Model_data_by_position %>%
        filter(task == "Base case") %>%
        reframe(percent = mean(percent, na.rm = TRUE),
                .by = c(rat_ID, rat_name, Genotype, Sex, 
                        task, detail, go, Response, Position)) %>%
        filter(Response == "Hit"))

Parametric_Check(Hit_position_AC_aov)

summary(Hit_position_AC_aov)

#### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "detail", "Position" # Main effects
), 
function(x) kruskal.test(reformulate(x, "percent"),
                         data = AC_Model_data_by_position %>%
                           filter(task == "Base case") %>%
                           reframe(percent = mean(percent, na.rm = TRUE),
                                   .by = c(rat_ID, rat_name, Genotype, Sex, 
                                           task, detail, go, Response, Position)) %>%
                           filter(Response == "Hit"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

#### Post-Hoc Dunn's Test ----
FSA::dunnTest(percent ~ interaction(Genotype, detail, Position),
              data = AC_Model_data_by_position %>%
                filter(task == "Base case") %>%
                reframe(percent = mean(percent, na.rm = TRUE),
                        .by = c(rat_ID, rat_name, Genotype, Sex, 
                                task, detail, go, Response, Position)) %>%
                filter(Response == "Hit"),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 3)[,1],
         detail1 = str_split_fixed(Comp1, '\\.', 3)[,2],
         position1 = str_split_fixed(Comp1, '\\.', 3)[,3],
         geno2 = str_split_fixed(Comp2, '\\.', 3)[,1],
         detail2 = str_split_fixed(Comp2, '\\.', 3)[,2],
         position2 = str_split_fixed(Comp1, '\\.', 3)[,3]) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(geno1 == geno2) %>%
  filter(position1 == position2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2) %>%
  print(n = 20)

## Reaction time -----
### Overall -----
Rxn_AC_aov = 
  aov(reaction ~ detail * Genotype,
      data = AC_Model_data %>%
        filter(task == "Base case") %>%
        filter(Response != "Miss") %>%
        reframe(reaction = mean(Rxn, na.rm = TRUE),
                .by = c(rat_ID, rat_name, Genotype, Sex, 
                        task, detail, go, Response)) %>%
        filter(Response == "Hit"))

Parametric_Check(Rxn_AC_aov)

summary(Rxn_AC_aov)

#### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "detail" # Main effects
), 
function(x) kruskal.test(reformulate(x, "reaction"),
                         data = AC_Model_data %>%
                           filter(task == "Base case") %>%
                           filter(Response != "Miss") %>%
                           reframe(reaction = mean(Rxn, na.rm = TRUE),
                                   .by = c(rat_ID, rat_name, Genotype, Sex, 
                                           task, detail, go, Response)) %>%
                           filter(Response == "Hit"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

#### Post-Hoc Dunn's Test ----
FSA::dunnTest(reaction ~ interaction(Genotype, detail),
              data = AC_Model_data %>%
                filter(task == "Base case") %>%
                filter(Response != "Miss") %>%
                reframe(reaction = mean(Rxn, na.rm = TRUE),
                        .by = c(rat_ID, rat_name, Genotype, Sex, 
                                task, detail, go, Response)) %>%
                filter(Response == "Hit"),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         detail1 = str_split_fixed(Comp1, '\\.', 2)[,2],
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         detail2 = str_split_fixed(Comp2, '\\.', 2)[,2]) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(geno1 == geno2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2)

#### KO vs. WT ----
wilcox.test(reaction ~ Genotype,
            data = AC_Model_data %>%
              filter(task == "Base case") %>%
              filter(Response == "Hit") %>%
              filter(detail == "Baseline") %>%
              reframe(reaction = mean(Rxn, na.rm = TRUE),
                      .by = c(rat_ID, rat_name, Genotype, Sex, 
                              task, detail, go, Response)) ) 


## Power analysis ----
### hit ----
AC_Model_data %>%
  filter(task == "Base case") %>%
  filter(Response == "Hit") %>%
  reframe(percent = mean(percent, na.rm = TRUE),
          .by = c(rat_ID, rat_name, Genotype, Sex, 
                  task, detail, Response)) %>%
  filter(detail %in% c("Baseline", "CNO 3mg/kg")) %>%
  mutate(level = factor(detail, labels = c(1, 2)) %>% as.numeric()) %>%
  select(level, percent) %>%
  # cramers_v()
  pearsons_c()

AC_hit_cohensD = 
lsr::cohensD(AC_Model_data %>%
               filter(task == "Base case") %>%
               filter(Response == "Hit") %>%
               reframe(percent = mean(percent, na.rm = TRUE),
                       .by = c(rat_ID, rat_name, Genotype, Sex, 
                               task, detail, Response)) %>%
               # filter(detail == "Between Treatment") %>%
               filter(detail == "Baseline") %>%
               .$percent,
             AC_Model_data %>%
               filter(task == "Base case") %>%
               filter(Response == "Hit") %>%
               reframe(percent = mean(percent, na.rm = TRUE),
                       .by = c(rat_ID, rat_name, Genotype, Sex, 
                               task, detail, Response)) %>%
               filter(detail == "CNO 3mg/kg") %>%
               .$percent)

pwr::pwr.t.test(d = AC_hit_cohensD,
                sig.level = 0.05, power = 0.80, type = "paired", alternative = "two.sided")
