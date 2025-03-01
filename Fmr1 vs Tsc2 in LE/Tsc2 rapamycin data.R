# Load Data (optional) ----------------------------------------------------

source("Fmr1 vs Tsc2 in LE BBN tone data.R")


# Get Treatment Group info ------------------------------------------------

## Group 1 ----
Tsc2_rapamycin_treated_rats_group1 = # c(315, 316, 318, 319, 306, 311)
  # dynamically selected based on treatment dates
  core_data %>%
  filter(line == "Tsc2-LE") %>% # only Tsc2 rats are being treated so this speeds up
  filter(date > 20240701 & date < 20241009) %>% # only Tsc2 rats are being treated so this speeds up
  filter(detail %in% c("Rapamycin (6mg/kg)")) %>%
  .$rat_ID %>% # use rat_ID because its unique
  unique # de-duplicate


## Group 2 ----
Tsc2_rapamycin_treated_rats_group2 = # c(328, 330:332, 370, 381, 384:385)
  # dynamically selected based on treatment dates
  core_data %>%
  filter(line == "Tsc2-LE") %>% # only Tsc2 rats are being treated so this speeds up
  filter(date > 20241001) %>% # only Tsc2 rats are being treated so this speeds up
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)")) %>%
  .$rat_ID %>% # use rat_ID because its unique
  unique # de-duplicate

# get list of rats with treatment
Tsc2_rapamycin_treated_rats = 
  core_data %>%
  filter(line == "Tsc2-LE") %>% # only Tsc2 rats are being treated so this speeds up
  filter(detail %in% c("Vehicle (Tween 80)", "Rapamycin (6mg/kg)",
                       "Post Treatment", "3+w Post Treatment",
                       "Rapamycin 2 (6mg/kg)", "Post Treatment 2")) %>% # treatment or control condition, to select for rats with both add filter(all(conditions))
  .$rat_ID %>% # use rat_ID because its unique
  unique # de-duplicate

Tsc2_treatment_dates =
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%
  arrange(date) %>% # force date order to be earliest 1st
  reframe(date = head(date, 1) %>% ymd(), #converted to date (year-month-day) for future math
          n = n(), # get count
          .by = c(rat_ID, detail)) %>% # for each task and detail
  arrange(rat_ID, date)


# Rxn data ----------------------------------------------------------------
Rapa_rxn_data =
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Recheck", "None",
                       "Vehicle (Tween 80)", "Post Vehicle", 
                       "Rapamycin (6mg/kg)", "Post Treatment", "3+w Post Treatment",
                       "Rapamycin 2 (6mg/kg)", "Post Treatment 2")) %>%
  mutate(detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
         detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)", "Post Vehicle",
                                            "Rapamycin (6mg/kg)", "Post Treatment", "3+w Post Treatment",
                                            "Rapamycin 2 (6mg/kg)", "Post Treatment 2"),
                         labels = c("Baseline", "Vehicle", "Post Vehicle", "Rapamycin",
                                    "Recovery", "Permanent",
                                    "Rapamycin 2", "Post Treatment 2"))) %>%
  filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  unnest(reaction) %>% 
  reframe(Rxn = mean(Rxn, na.rm = TRUE) * 1000,
          .by = c(rat_ID, rat_name, sex, genotype, line, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`))

## 2 weeks prior to treatment ----
Rapa_rxn_data_limited =
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Recheck", "None",
                       "Vehicle (Tween 80)", "Post Vehicle", 
                       "Rapamycin (6mg/kg)", "Post Treatment", "3+w Post Treatment",
                       "Rapamycin 2 (6mg/kg)", "Post Treatment 2")) %>%
  mutate(detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
         detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)", "Post Vehicle",
                                            "Rapamycin (6mg/kg)", "Post Treatment", "3+w Post Treatment",
                                            "Rapamycin 2 (6mg/kg)", "Post Treatment 2"),
                         labels = c("Baseline", "Vehicle", "Post Vehicle", "Rapamycin",
                                    "Recovery", "Permanent",
                                    "Rapamycin 2", "Post Treatment 2"))) %>%
  filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  # filter to only 2 weeks prior to treatment
  left_join(Tsc2_treatment_dates %>% 
              filter(detail == "Vehicle (Tween 80)") %>% 
              reframe(start_date = (date - 14) %>% str_remove_all(pattern = "-"), .by = rat_ID),
            by = join_by(rat_ID)) %>%
  rowwise() %>%
  filter(date > start_date) %>%
  # calculate reaction time average
  unnest(reaction) %>% 
  reframe(Rxn = mean(Rxn, na.rm = TRUE) * 1000,
          .by = c(rat_ID, rat_name, sex, genotype, line, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`))

