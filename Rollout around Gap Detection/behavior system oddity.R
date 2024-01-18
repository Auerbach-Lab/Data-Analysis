
# Get rats running around GD roll out -------------------------------------

# roll out date = 1/10/2023

rats_at_rollout = 
  intersect(filter(run_archive, date == 20230109) %>% reframe(rat_name = unique(rat_name), rat_ID = unique(rat_ID)),
            filter(run_archive, date == 20230111) %>% reframe(rat_name = unique(rat_name), rat_ID = unique(rat_ID)))



# Prep data ---------------------------------------------------------------

rollout_data_nested = 
  run_archive %>%
  # select rats active at before and after roll out
  filter(rat_ID %in% rats_at_rollout$rat_ID) %>%
  # select only dates around roll out 
  filter(date > 20221101 & date < 20230301 & date != 20230103 & date != 20230104) %>%
  mutate(rollout = if_else(date < 20230110, "Pre", "Post"),
         rollout = factor(rollout, levels = c("Pre", "Post"))) %>%
  unnest_wider(assignment) %>% 
  arrange(rat_name, date)


rollout_data_filtered = 
  rollout_data_nested %>%
  # drop Oddball because it can't be filtered by Frequency since on rolling frequency
  filter(experiment != "Oddball") %>%
  group_by(rat_ID) %>%
  do(as_tibble(.) %>% 
       mutate(Frequency = str_extract(file_name, pattern = ".*?(?=_)")) %>%
       # why are you not actually filtering?????
       rowwise %>%
       filter(Frequency == filter(., date == "20230110") %>% select(Frequency) %>% unique())) %>% 
  mutate(Frequency = factor(Frequency, levels = c("4kHz", "8kHz", "16kHz", "32kHz", "BBN"))) %>%
  # select(1, 4, 5, 7:10, Frequency) %>% 
  unnest_wider(stats)


rollout_data1 = 
  rollout_data_filtered %>%
  select(-dprime) %>%
  unnest(reaction, names_sep = ".")

rollout_data2 =
  rollout_data_filtered %>%
  select(date, time, rat_name, rat_ID, file_name, dprime) %>%
  unnest(dprime, names_sep = ".", keep_empty = TRUE)

rollout_data =
  left_join(rollout_data1,
            rollout_data2,
            by = join_by(date, time, rat_name, rat_ID, file_name,
                         # columns that are shared but have unique names
                         `reaction.Freq (kHz)` == dprime.Freq, `reaction.Inten (dB)` == dprime.dB, `reaction.Dur (ms)` == dprime.Dur)) %>% 
  rename(kHz = `reaction.Freq (kHz)`,
         dB = `reaction.Inten (dB)`,
         dur = `reaction.Dur (ms)`,
         reaction = `reaction.Rxn`,
         dprime = `dprime.dprime`)



# Basic Stats Graphs ------------------------------------------------------

General_Stats =
  rollout_data_filtered %>%
    select(date, rat_ID, rat_name, trial_count, hit_percent, FA_percent, rollout) %>%
    mutate(hit_percent = hit_percent * 100,
           FA_percent = FA_percent * 100) %>%
    unique() %>%
    ungroup %>%
    reframe(trial_count = mean(trial_count, na.rm = TRUE),
            hit_percent = mean(hit_percent, na.rm = TRUE),
            FA_percent = mean(FA_percent, na.rm = TRUE),
            .by = c(rat_ID, rat_name, rollout)) %>%
  pivot_longer(
    cols = c(trial_count, hit_percent, FA_percent), 
    names_to = "Stat",
    values_to = "value"
  ) %>%
  mutate(Stat = factor(Stat, levels = c("trial_count", "hit_percent", "FA_percent"))) %>%
  ggplot(aes(x = rollout, y = value, color = rollout, group = rat_ID)) +
    geom_boxplot(aes(group = rollout, fill = Stat), show.legend = FALSE, linewidth = 1) +
    # stat_summary(fun.data = boxplot_quantiles, geom = "boxplot", na.rm = TRUE, linewidth = 1) +
    # stat_summary(fun = mean,
    #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
    #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
    #              geom = "errorbar", width = 0, position = position_dodge(1)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
                 geom = "line", color = "black") +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 geom = "point") +
    scale_fill_manual(values = c("lightgreen", "darkgrey", "tan"))+
    facet_wrap(~ Stat, scales = "free") +
    theme_ipsum_es()

print(General_Stats)

ggsave(filename = "General_Stats.jpg",
       path = save_folder,
       plot = General_Stats,
       width = 12, height = 6, units = "in", dpi = 300)


# Reaction Curves Pre vs. Post --------------------------------------------
# Not including Oddball

Rxn_all_frequencies = 
  rollout_data %>%
    filter(dur == 50) %>%
    ggplot(aes(x = dB, y = reaction * 1000, color = rollout, group = rollout)) +
    # stat_summary(fun = mean,
    #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
    #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
    #              geom = "errorbar", width = 0, position = position_dodge(1)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
                 geom = "line") +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 geom = "point", size = 2) +
    scale_x_continuous(breaks = seq(0, 100, by = 20)) +
    facet_wrap( ~ Frequency, ncol = 2) +
    theme_ipsum_es()

print(Rxn_all_frequencies)

ggsave(filename = "Rxn_curve_all_frequencies.jpg",
       path = save_folder,
       plot = Rxn_all_frequencies,
       width = 12, height = 6, units = "in", dpi = 300)


Rxn_BNN_all_durations = 
  rollout_data %>%
  filter(Frequency == "BBN") %>%
  ggplot(aes(x = dB, y = reaction * 1000, color = rollout, group = rollout)) +
  # stat_summary(fun = mean,
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 0, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", size = 2) +
  
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  facet_wrap( ~ dur, ncol = 3) +
  theme_ipsum_es()

print(Rxn_BNN_all_durations)

ggsave(filename = "Rxn_curve_all_durations.jpg",
       path = save_folder,
       plot = Rxn_BNN_all_durations,
       width = 12, height = 6, units = "in", dpi = 300)


# Individual Graphs
ggplot(data = rollout_data,
       aes(x = dB, y = reaction, color = rollout, group = interaction(dur, rollout))) +
  geom_line(aes(color = rollout, group = interaction(date, dur)),
            alpha = 0.5) +
  geom_smooth(se = FALSE, linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  xlim(20, 60) +
  facet_wrap(~ rat_ID, scales = "free") +
  theme_ipsum_es()


# Select Reaction times Pre vs. Post --------------------------------------

Rxn_boxplot_all_durations = 
  rollout_data %>%
    filter(dB %in% c(30, 60)) %>%
    ggplot(aes(x = interaction(rollout, dB), y = reaction * 1000, color = rollout, group = interaction(dB, rat_ID))) +
    geom_boxplot(aes(fill = dB, group = interaction(rollout, dB)), na.rm = TRUE, linewidth = 1) +
    # stat_summary(fun = mean,
    #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
    #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
    #              geom = "errorbar", width = 0, position = position_dodge(1)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
                 geom = "line", color = "goldenrod") +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 geom = "point") +
    facet_wrap( ~ dur) +
    theme_ipsum_es()

print(Rxn_boxplot_all_durations)

ggsave(filename = "Rxn_boxplot_all_durations.jpg",
       path = save_folder,
       plot = Rxn_boxplot_all_durations,
       width = 12, height = 6, units = "in", dpi = 300)

# Thresholds by Duration, Pre vs. Post ------------------------------------

TH_daily_boxplot = 
  rollout_data_filtered %>%
    select(date, rat_ID, rat_name, file_name, experiment, phase, task, detail, rollout, threshold) %>%
    unique %>%
    unnest(threshold) %>% 
    filter(! is.na(TH)) %>% 
    ggplot(aes(x = rollout, y = TH, color = rat_name, group = rat_ID)) +
      stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
                   geom = "line", color = "black") +
      stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                   geom = "point", size = 4) +
      facet_wrap(~ Dur) +
      theme_ipsum_es()

print(TH_daily_boxplot)

ggsave(filename = "TH_daily_boxplot.jpg",
       path = save_folder,
       plot = TH_daily_boxplot,
       width = 12, height = 6, units = "in", dpi = 300)


# dprime curve graphs -----------------------------------------------------

dprime_all_frequencies = 
  rollout_data %>%
  filter(dur == 50) %>%
  ggplot(aes(x = dB, y = dprime, color = rollout, group = rollout)) +
  geom_hline(yintercept  = 1.5, linetype = "longdash", color = "darkblue") +
  # stat_summary(fun = mean,
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 0, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", size = 2) +
  facet_wrap( ~ Frequency, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme_ipsum_es()

print(dprime_all_frequencies)

ggsave(filename = "dprime_curve_all_frequencies_50ms_duration.jpg",
       path = save_folder,
       plot = dprime_all_frequencies,
       width = 12, height = 6, units = "in", dpi = 300)


dprime_BNN_all_durations = 
  rollout_data %>%
  filter(Frequency == "BBN") %>%
  ggplot(aes(x = dB, y = dprime, color = rollout, group = rollout)) +
  geom_hline(yintercept  = 1.5, linetype = "longdash", color = "darkblue") +
  # stat_summary(fun = mean,
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 0, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", size = 2) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  facet_wrap( ~ dur, ncol = 3) +
  theme_ipsum_es()

print(dprime_BNN_all_durations)

ggsave(filename = "dprime_curve_all_durations.jpg",
       path = save_folder,
       plot = dprime_BNN_all_durations,
       width = 12, height = 6, units = "in", dpi = 300)


# Individual Graphs
ggplot(data = rollout_data,
       aes(x = dB, y = dprime, color = rollout, group = interaction(dur, rollout))) +
  geom_line(aes(color = rollout, group = interaction(date, dur)),
            alpha = 0.5) +
  geom_smooth(se = FALSE, linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  xlim(20, 40) +
  facet_wrap(~ rat_ID, scales = "free") +
  theme_ipsum_es()


  