# Receiver Operating Characteristic (ROC) Curve
# X = FA%   (0-1)
# y = hit%  (0-1)
# Chance line: Hit% = FA%
# Neutral Bias (1/2 diagonal): Hit% + FA% = 1
# # below = conservative, above = liberal

# TODO: graph for each individual before in:
#              - different backgrounds
#              - near threshold vs 30-90dB


# use stats_table, stats_table_by_BG & stats_table_detail



# Custom graph bits ------------------------------------------------------------

bias_theme = list(
  scale_fill_manual(values = c("baseline" = "black", "HL" = "blue", "recovery" = "yellow", "post-HL" = "green"),
                    labels = c("baseline" = "Baseline", "recovery" = "Recovery", "post-HL" = "Post sound\nexposure\n(3+ weeks)")),
  scale_x_continuous(expand = c(0, 0.01)),
  scale_y_continuous(expand = c(0, 0.01)),
  coord_cartesian(xlim = c(0, 0.35),
                  ylim = c(0.65, 1)),
  theme_classic(),
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ))

divider_line = list(geom_segment(aes(x = 0, y = 1, xend = 0.5, yend = 0.5), color = "darkgreen", linewidth = 1))

chance_line = list(geom_abline(intercept = 0, slope = 1, color = "darkred", linewidth = 1))

# Starter Graph -----------------------------------------------------------
# Includes everything (background, mixed presentation, CNO treatment, etc)
stats_average = stats_average %>%
  reframe(hit_percent = mean(hit_percent, na.rm = TRUE),
          FA_percent = mean(FA_percent, na.rm = TRUE),
          .by = c(HL_state, stim_type))
    
stats_table %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state != "HL") %>%
  ggplot(aes(x = FA_percent, y = hit_percent, color = stim_type, fill = HL_state)) +
  divider_line +
  geom_path(data = stats_average %>% filter(HL_state != "HL"), 
            aes(group = interaction(stim_type)), alpha = 0.8, linewidth = 1.5) + 
  geom_path(aes(group = interaction(rat_ID, stim_type, Sex))) +
  geom_point(shape = 21, size = 2) +
  geom_point(data = stats_average %>% filter(HL_state != "HL"),
             size = 6, alpha = 0.7, stroke = 2, shape = 21, show.legend = FALSE) +
  labs(title = "Individual Bias",
       x = "False Alarm %", y = "Hit %",
       color = "Sound Type", fill = "Time point"
       ) +
  bias_theme

# Bias with Sex Graph -----------------------------------------------------------
stats_average_sex = stats_table %>%
  reframe(hit_percent = mean(hit_percent, na.rm = TRUE),
          FA_percent = mean(FA_percent, na.rm = TRUE),
          .by = c(HL_state, stim_type, Sex))


stats_table %>%
  # filter(stim_type == "tone") %>%
  mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state != "HL") %>%
  ggplot(aes(x = FA_percent, y = hit_percent, shape = Sex, color = stim_type, fill = HL_state)) +
  geom_path(data = stats_average_sex %>% filter(HL_state != "HL"), 
            aes(group = interaction(stim_type, Sex)), alpha = 0.8, linewidth = 1.5) + 
  geom_segment(aes(x = 0, y = 1, xend = 0.5, yend = 0.5), color = "darkgreen", linewidth = 1) +
  geom_path(aes(group = interaction(rat_ID, stim_type, Sex))) +
  geom_point(size = 2) +
  stat_summary(data = stats_average_sex %>% filter(HL_state != "HL"), 
               show.legend = FALSE, size = 6, alpha = 0.7, stroke = 2,
               fun = "mean", na.rm = TRUE, geom = "point") +
  geom_text(label="Convervative", color = "black", size = 10,
            x = 0.1, y = 0.68) +
  geom_text(label="Libral", color = "black", angle = 90, size = 10,
            x = 0.33, y = 0.85) +
  scale_shape_manual(values = c("Male" = 21, "Female" = 24)) +
  scale_fill_manual(values = c("baseline" = "black", "HL" = "blue", "recovery" = "yellow", "post-HL" = "green"),
                    labels = c("baseline" = "Baseline", "recovery" = "Recovery", "post-HL" = "Post sound\nexposure\n(3+ weeks)")) +
  labs(title = "Bias by Sex",
       x = "False Alarm %", y = "Hit %",
       color = "Sound Type", fill = "Time point", shape = "Sex",
  ) +
  scale_x_continuous(expand = c(0, 0.01)) +
  scale_y_continuous(expand = c(0, 0.01)) +
  coord_cartesian(xlim = c(0, 0.35),
                  ylim = c(0.65, 1)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


# Data Gathering and filtering --------------------------------------------

## Get only appropriate trails for each rat --
# This includes CNO treatment days, BG noise, Mixed BBN and tone days
core_runs_list =
  core_data %>% select(rat_name, rat_ID, UUID, Sex, HL_state, BG_type, BG_Intensity, date, file_name, phase, task, detail) %>%
  mutate(file_dB_range = str_extract(file_name, pattern = "[:digit:]+?-*[:digit:]+?dB|MIX.*?dB") %>% str_remove(pattern = "dB$"))

trial_table = trial_archive %>%
  filter(UUID %in% core_runs_list$UUID) %>%
  left_join(core_runs_list, by = join_by(UUID)) %>%
  select(-UUID)

# BBN ---------------------------------------------------------------------

## Data ----
bias_table_BBN = trial_table %>%
  #only days with BBN presentation
  filter(`Freq (kHz)` == 0) %>%
  # drop Mixed days because they only have a 300ms No Go even though there are 3 durations of Go stimuli
  filter(detail != "Mixed") %>%
  summarise(hit_percent = sum(Response == "Hit") / sum(Response %in% c("Hit", "Miss")),
            FA_percent = sum(Response == "FA") / sum(Response %in% c("CR", "FA")),
            hit = sum(Response == "Hit"),
            Go = sum(Response %in% c("Hit", "Miss")),
            FA = sum(Response == "FA"),
            NG =  sum(Response %in% c("CR", "FA")),
            .by = c(rat_name, rat_ID, Sex, HL_state, BG_type, BG_Intensity, `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`))

### BBN no BG ------
bias_table_BBN_filtered =
  bias_table_BBN  %>%
    # only days with no Background
    filter(BG_type == "None") %>%
    filter(HL_state != "HL") %>%
    # keep only FAs and easiest intensity that is shared which is 80 NOT 90
    filter(`Inten (dB)` %in% c(20, 30)) %>%
    select(-FA_percent, -FA, -NG) %>%
    left_join(bias_table_BBN  %>%
                filter(HL_state != "HL") %>%
                # keep only FAs and easiest intensity that is shared which is 80 NOT 90
                filter(`Inten (dB)` %in% c(-100)) %>% select(rat_name:`Freq (kHz)`, `Dur (ms)`, FA_percent, FA, NG),
              by = join_by(rat_name, rat_ID, Sex, HL_state, BG_type, BG_Intensity, `Freq (kHz)`, `Dur (ms)`)) %>%
    mutate(`Dur (ms)` = as.factor(`Dur (ms)`))

bias_table_BBN_average =  bias_table_BBN_filtered %>%
    summarise(hit_percent = mean(hit_percent, na.rm = TRUE),
              FA_percent = mean(FA_percent, na.rm = TRUE),
              .by = c(HL_state, `Dur (ms)`, `Inten (dB)`)) %>%
  mutate(`Dur (ms)` = as.factor(`Dur (ms)`)) %>%
  arrange(HL_state)

### BBN BG ------
# Doesn't exist only did tones with and without background

## Graphs -----
    
### BBN no BG bias graph ------ 
bias_table_BBN_filtered %>%
  ggplot(aes(x = FA_percent, y = hit_percent, fill = HL_state, shape = `Dur (ms)`, group = `Dur (ms)`)) +
    chance_line +
    divider_line +
    geom_path(data = bias_table_BBN_average, alpha = 0.8, linewidth = 1.5, show.legend = FALSE) +
    geom_path(aes(group = interaction(rat_ID, `Dur (ms)`))) +
    geom_point(size = 2, stroke = 1) +
    geom_point(data = bias_table_BBN_average, size = 5, stroke = 2, alpha = 0.6, show.legend = FALSE) +
    labs(title = "Individual bias near threshold",
         x = "False Alarm %", y = "Hit %",
         color = "Sound Type", fill = "Time point"
    ) +
    bias_theme +
    facet_wrap(~ `Inten (dB)`) +
    # labels
    geom_text(label = "Convervative", color = "black", size = 10,
              x = 0.06, y = 0.12) +
    geom_text(label = "Libral", color = "black", angle = 90, size = 10,
              x = 0.45, y = 0.85) +
    geom_text(label = "Guessing/\nChance", color = "darkred", size = 5.5,
              x = 0.35, y = 0.25) +
    scale_shape_manual(values = c("300" = 21, "100" = 24, "50" = 23)) +
    coord_cartesian(xlim = c(0, 0.5),
                    ylim = c(0.1, 1))


# Tones -------------------------------------------------------------------


### Single frequency per file ------
# This table is challenging to graph since FA only happen on -100
bias_table_quiet_frequency = trial_table %>%
  # only days with no Background
  filter(BG_type == "None") %>%
  summarise(hit_percent = sum(Response == "Hit") / sum(Response %in% c("Hit", "Miss")),
            FA_percent = sum(Response == "FA") / sum(Response %in% c("CR", "FA")),
            hit = sum(Response == "Hit"),
            Go = sum(Response %in% c("Hit", "Miss")),
            FA = sum(Response == "FA"),
            NG =  sum(Response %in% c("CR", "FA")),
            .by = c(rat_name, rat_ID, Sex, HL_state, BG_type, BG_Intensity, `Freq (kHz)`, file_dB_range, `Dur (ms)`))

# Alternatively could sort into Rxn vs TH days task days
easy_range = c("25-85", "30-90", "60-90", "50-80", "30-60")
hard_range = c("15-75", "10-70", "20-50", "10-40", "0-60", "25-55", "30-60", "MIX", "MIX5step", "15-45", "10-80")
not_used_ranges = c("20-80", "30-60", "40-70", "45-75", "30")

bias_table_quiet_binned = 
  trial_table %>%
  # only days with no Background
  filter(BG_type == "None") %>%
  #only days with single frequency presentation
  filter(detail == "Alone") %>%
  summarise(hit_percent_easy = sum(Response == "Hit" & file_dB_range %in% easy_range) / 
              sum(Response %in% c("Hit", "Miss") & file_dB_range %in% easy_range),
            FA_percent_easy = sum(Response == "FA" & file_dB_range %in% easy_range) / 
              sum(Response %in% c("CR", "FA") & file_dB_range %in% easy_range),
            hit_percent_hard = sum(Response == "Hit" & file_dB_range %in% hard_range) / 
              sum(Response %in% c("Hit", "Miss") & file_dB_range %in% hard_range),
            FA_percent_hard = sum(Response == "FA" & file_dB_range %in% hard_range) / 
              sum(Response %in% c("CR", "FA") & file_dB_range %in% hard_range),
            hit_easy = sum(Response == "Hit" & file_dB_range %in% easy_range),
            Go_easy = sum(Response %in% c("Hit", "Miss") & file_dB_range %in% easy_range),
            FA_easy = sum(Response == "FA" & file_dB_range %in% easy_range),
            NG_easy = sum(Response %in% c("CR", "FA") & file_dB_range %in% easy_range),
            hit_hard = sum(Response == "Hit" & file_dB_range %in% hard_range),  
            Go_hard = sum(Response %in% c("Hit", "Miss") & file_dB_range %in% hard_range),
            FA_hard = sum(Response == "FA" & file_dB_range %in% hard_range),
            NG_hard = sum(Response %in% c("CR", "FA") & file_dB_range %in% hard_range),
            .by = c(rat_name, rat_ID, Sex, HL_state, BG_type, BG_Intensity, `Freq (kHz)`, `Dur (ms)`)) %>% View


# Get Averages
reframe(hit_percent = mean(hit_percent, na.rm = TRUE),
        FA_percent = mean(FA_percent, na.rm = TRUE),
        .by = c(rat_ID, rat_name, Sex, HL_state, stim_type, Freq)) %>% View

mutate(HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state != "HL") %>%
  ggplot(aes(x = FA_percent, y = hit_percent, color = stim_type, fill = HL_state)) +
  # geom_abline(intercept = 0, slope = 1, color = "darkgreen", linewidth = 1) +
  geom_segment(aes(x = 0, y = 1, xend = 0.5, yend = 0.5), color = "darkgreen", linewidth = 1) +
  geom_line(aes(group = interaction(rat_ID, stim_type))) +
  geom_point(shape = 21, size = 2)




# Individual change in hit and FA across HL_state -------------------------
