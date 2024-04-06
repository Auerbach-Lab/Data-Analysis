# All genotypes graph ----------
Rxn_table %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(Duration %in% c(50, 100, 300)) %>%
  # rename(Intensity = `Inten (dB)`) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  # filter(rat_ID < 314) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == single_Frequency) %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(group), shape = line,
             color = genotype, group = interaction(Duration, group, line, genotype))) +
  ## Overall average lines
  # stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
  #              fun = function(x) mean(x, na.rm = TRUE),
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 1.5, linetype = "solid", linewidth = 1.5, alpha = 0.5,
  #              position = position_dodge(4)) +
  # stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
  #              fun = function(x) mean(x, na.rm = TRUE),
  #              geom = "line", linetype = "solid", linewidth = 1.5, alpha = 0.5,
  #              position = position_dodge(4)) +
## Lines for each group
stat_summary(fun = function(x) mean(x, na.rm = TRUE),
             fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
             fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
             geom = "errorbar", width = 1.5, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line", position = position_dodge(1)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "", shape = "ASD model",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
  scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 


# LE Data Loading ---------------------------------------------------------

dataset_LE = run_archive %>%
  # Omit Invalid runs
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a tibble
  suppressWarnings(
    filter(str_detect(as.vector(warnings_list), pattern = "wrong delay window", negate = TRUE)) 
  ) %>%
  # Only include LE rats
  filter(rat_ID %in% c(265, 266, 267, 268)) %>%
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>%
  # Only keep relevant Experiments
  filter(experiment %in% c("TTS")) %>%
  filter(! phase %in% c("Octave")) %>%
  filter(! task %in% c("Duration Testing")) %>%
  filter(! task %in% c("Training", "Reset")) %>%
  unnest_wider(stats) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # record date of hearing loss
  left_join(select(rat_archive, Rat_ID, Genotype, Sex), by = c("rat_ID" = "Rat_ID")) %>%
  # effectively hardcoded
  mutate(line = Genotype,
         genotype = "WT") %>%
  # Calculate time since HL
  mutate(BG = if_else(str_detect(file_name, pattern = "_BG_"), str_extract(file_name, pattern = "BG_.*$"), "None")) %>%
  filter(BG == "None") %>%
  rename(sex = Sex)

TH_table_LE =
  dataset_LE %>%
  # re-nest by Frequency
  unnest(dprime) %>%
  filter(Freq == "0") %>%
  nest(dprime = c(rat_name, date, Freq, Dur, dB, dprime), 
       .by = c(rat_ID, rat_name, Freq, Dur, line, genotype, sex, detail)) %>%
  # calculate TH
  rowwise() %>%
  mutate(TH = Calculate_TH(dprime)) %>%
  select(-dprime) %>%
  rename(Frequency = Freq, Duration = Dur)

Rxn_table_LE =
  dataset_LE %>%
  # Get Reaction times:
  unnest(reaction) %>%
  rename(Frequency = `Freq (kHz)`, Duration = `Dur (ms)`, Intensity = `Inten (dB)`) %>%
  filter(Frequency == "0") %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, Frequency, Duration, Intensity, sex, genotype, line, detail) %>%
  # Get Averages
  do(transmute(., Rxn = mean(Rxn, na.rm = TRUE) * 1000)) %>%
  unique()


# WT & LE Comparisons -----------------------------------------------------


## TH Averages -------------------------------------------------------------
TH_table %>%
  filter(Frequency == 0) %>%
  filter(genotype == "WT") %>%
  bind_rows(TH_table_LE) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  group_by(line, genotype, detail, Frequency , Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")

# TH Graph ----------------------------------------------------------------

TH_table %>%
  filter(Frequency == 0) %>%
  filter(genotype == "WT") %>%
  bind_rows(TH_table_LE) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(! all_of(group == "Group 2" & detail == "Alone" & line == "Tsc2-LE")) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
           factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH, shape = line,
             fill = detail, color = line, group = interaction(detail, line, genotype))) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
               position = position_dodge(1), vjust = 2, size = 3) +
  # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
  scale_fill_manual(values = c("Alone" = "darkgrey", "Recheck" = "goldenrod")) +
  scale_color_manual(values = c("LE" = "black", "Tsc2-LE" = "deepskyblue", "Fmr1-LE" = "lightcoral")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3"),
       fill = "Genotype") +
  facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

# WT & LE graph ----------

Rxn_table %>%
  filter(Frequency == 0) %>%
  filter(genotype == "WT") %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(Duration %in% c(50, 100, 300)) %>%
  bind_rows(Rxn_table_LE) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  # filter(rat_ID < 314) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(genotype == "WT") %>%
  filter(! all_of(group == "Group 2" & detail == "Alone" & line == "Tsc2-LE" & Duration %in% c(50, 300))) %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(group), shape = line, fill = line,
             color = genotype, group = interaction(Duration, group, line, genotype))) +
  ## Overall average lines
  # stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
  #              fun = function(x) mean(x, na.rm = TRUE),
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 1.5, linetype = "solid", linewidth = 1.5, alpha = 0.5,
  #              position = position_dodge(4)) +
  # stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
  #              fun = function(x) mean(x, na.rm = TRUE),
  #              geom = "line", linetype = "solid", linewidth = 1.5, alpha = 0.5,
  #              position = position_dodge(4)) +
## Lines for each group
stat_summary(fun = function(x) mean(x, na.rm = TRUE),
             fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
             fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
             geom = "errorbar", width = 1.5, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line", position = position_dodge(1)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "", shape = "Line", fill = "Line",
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
  scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_shape_manual(values = c("LE" = 21, "Fmr1-LE" = 22, "Tsc2-LE" = 24)) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 
