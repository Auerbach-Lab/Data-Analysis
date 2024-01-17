
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
  arrange(rat_name, date) %>%
  unnest_wider(stats)


rollout_data1 = 
  rollout_data_nested %>%
  select(-dprime) %>%
  unnest(reaction, names_sep = ".")

rollout_data2 =
  rollout_data_nested %>%
  select(-reaction) %>%
  unnest(dprime, names_sep = ".", keep_empty = TRUE)

rollout_data =
  left_join(rollout_data1,
            rollout_data2 %>%
              # remove Oddball which is missing the d' data
              filter(experiment %in% c("Oddball", "Oddball Training")),
            by = join_by(date, time, box, rat_name, rat_ID, weight, 
                         file_name, assigned_file_name, experiment, phase, task, detail, comment, 
                         summary, stim_type, analysis_type, 
                         trial_count, hits, misses, CRs, FAs, hit_percent, FA_percent, mean_attempts_per_trial, 
                         threshold, FA_detailed, block_size, complete_block_count,
                         invalid, comments, warnings_list, omit_list, UUID, scientist, weightProblem, rxnProblem,
                         rollout,
                         # columns that are shared but have unique names
                         `reaction.Freq (kHz)` == dprime.Freq, `reaction.Inten (dB)` == dprime.dB, `reaction.Dur (ms)` == dprime.Dur)) %>%
  rename(kHz = `reaction.Freq (kHz)`,
         dB = `reaction.Inten (dB)`,
         dur = `reaction.Dur (ms)`,
         reaction = `reaction.Rxn`,
         dprime = `dprime.dprime`)


# Look at files during the rollout period ---------------------------------

rollout_data_filtered = 
  rollout_data %>%
    filter(experiment != "Oddball") %>%
    group_by(rat_ID) %>%
    do(as_tibble(.) %>% 
         # .[c(1, 4, 5, 7:10)] %>%
         mutate(Freq = str_extract(file_name, pattern = ".*?(?=_)")) %>%
         # why are you not actually filtering?????
         rowwise %>%
         filter(phase == (filter(., date == "20230110") %>% select(phase) %>% unique()) & 
                  Freq == filter(., date == "20230110") %>% select(Freq) %>% unique()))



# Reaction Curves Pre vs. Post --------------------------------------------
# Not including Oddball
ggplot(data = rollout_data_filtered,
       aes(x = dB, y = reaction, color = rollout, group = interaction(dur, rollout))) +
  geom_line(aes(color = rollout, group = interaction(date, dur)),
            alpha = 0.5) +
  geom_smooth(se = FALSE, linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  xlim(20, 60) +
  facet_wrap(~ rat_ID, scales = "free") +
  theme_ipsum_es()


# Select Reaction times Pre vs. Post --------------------------------------

rollout_data_filtered %>%
  filter(dB %in% c(30, 35, 40, 50)) %>%
  ggplot(aes(x = rollout, y = reaction, color = rat_name, group = interaction(dB, rat_ID))) +
  # stat_summary(fun = mean,
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 0, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", size = 2) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line") +
  facet_wrap(~ dur, scales = "free") +
  theme_ipsum_es()
  


# Thresholds by Duration, Pre vs. Post ------------------------------------

rollout_data_filtered %>%
  select(date, rat_ID, rat_name, file_name, experiment, phase, task, detail, rollout, threshold) %>%
  unique %>%
  unnest(threshold) %>% 
  filter(! is.na(TH)) %>% 
  ggplot(aes(x = rollout, y = TH, color = rat_name, group = rat_ID)) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                 geom = "point", size = 2) +
    stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
                 geom = "line") +
    facet_wrap(~ Dur, scales = "free") +
    theme_ipsum_es()
  
  