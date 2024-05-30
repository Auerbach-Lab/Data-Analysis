# Load Data (optional) ----------------------------------------------------

# source("Fmr1 vs Tsc2 in LE BBN tone data.R")




# Get Treatment Group info ------------------------------------------------

# get list of rats with treatment
Tsc2_rapamycin_treated_rats = 
  core_data %>%
  filter(line == "Tsc2-LE") %>% # only Tsc2 rats are being treated so this speeds up
  filter(detail %in% c("Vehicle (Tween 80)", "Rapamycin")) %>% # treatment or control condition, to select for rats with both add filter(all(conditions))
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


# Piloting data checks ----------------------------------------------------

TH_table %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(Duration == 50) %>%
  filter(detail %in% c("Recheck", "Vehicle (Tween 80)", "Rapamycin")) %>%
  reframe(TH = mean(TH, na.rm = TRUE),
          .by = c(genotype, detail))

vehicle_check_rxn =
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Recheck", "Vehicle (Tween 80)", "Rapamycin")) %>%
  filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  left_join(Tsc2_treatment_dates %>% 
              filter(detail == "Vehicle (Tween 80)") %>% 
              reframe(start_date = (date - 14) %>% str_remove_all(pattern = "-"), .by = rat_ID),
            by = join_by(rat_ID)) %>%
  rowwise() %>%
  filter(date > start_date) %>%
  unnest(reaction) %>% 
  reframe(Rxn = mean(Rxn, na.rm = TRUE) * 1000,
          .by = c(rat_ID, rat_name, sex, genotype, line, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`))

## Rxn Graph ====
vehicle_check_rxn %>%
  # filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  filter(`Inten (dB)` > 15) %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn, linetype = as.factor(detail),
             color = genotype, group = interaction(detail, genotype))) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 1, position = position_dodge(0.5)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(0.5)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", position = position_dodge(0.5)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
  # scale_linetype_manual(values = c("Recheck" = "solid", "Vehicle (Tween 80)" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 
