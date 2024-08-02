# Load Necessary Datasets -------------------------------------------------

rat_genotypes = fread(glue("{ABR_data_folder}/Genotypes.csv"))


# Load CSV file generated from the abr_processor script (credit: Walker Gauthier)
ABR_single_csv = fread(glue("{ABR_data_folder}/ABR_data.csv"))

ABR_data = left_join(ABR_single_csv %>%
                       as_tibble() %>%
                       # drop Wave not picked 
                       filter(Wave_amp != 0 & Wave_lat != 0) %>%
                       mutate(rat_ID = str_extract(file, pattern = "00[:digit:]+") %>% as.numeric()),
                     rat_genotypes,
                     by = join_by(rat_ID))



# Write our data ----------------------------------------------------------
fwrite(ABR_data, glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Walker/Fmr1_SD_ABR_data", str_remove_all(Sys.Date(), "-"),".csv"), row.names = FALSE)


# Prep data ---------------------------------------------------------------

# Ensure one entry for each rat, ear, frequency & intensity
# Each ear is a separate n
ABR_data_ear = 
  ABR_data %>%
    summarise(Wave_amp = mean(Wave_amp, na.rm = TRUE),
              Wave_lat = mean(Wave_lat, na.rm = TRUE),
              RMS = mean(current_rms, na.rm = TRUE),
              .by = c(rat_ID, genotype, current_ear, current_freq, current_wave, current_inten))

ABR_data_rat = 
  ABR_data %>%
  summarise(Wave_amp = mean(Wave_amp, na.rm = TRUE),
            Wave_lat = mean(Wave_lat, na.rm = TRUE),
            RMS = mean(current_rms, na.rm = TRUE),
            .by = c(rat_ID, genotype, current_freq, current_wave, current_inten))

## Sample counts =====
ABR_data_rat %>% 
  select(rat_ID, genotype) %>% 
  unique %>% 
  reframe(n = n(), .by = c(genotype))

# Ben Poster Plots --------------------------------------------------------

## Amplitude =====
Ben_Amp_plot =
ggplot(data = filter(ABR_data_rat, current_inten == 80),
       aes(x = current_wave, y = Wave_amp, fill = genotype, group = interaction(current_wave, genotype, current_wave))) +
  geom_boxplot(linewidth = 1) +
  # geom_point(show.legend = FALSE, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "Wave",
       y = "Amplitude @ 80dB",
       fill = "Genotype") +
  facet_wrap(~ current_freq, scales = "free_y") +
  theme_classic(base_size = 20) +
  theme(
    # legend.position=c(.8,.8)
    legend.position=c(.8,.2)
  )

print(Ben_Amp_plot)

ggsave(file="Fmr1_SD_ABR_amp.svg", path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben",
       plot = Ben_Amp_plot, 
       width = 10, height = 8, dpi = 600)

ggsave(file="Fmr1_SD_ABR_amp.jpg", path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben",
       plot = Ben_Amp_plot, 
       width = 10, height = 8, dpi = 600)

## Latency =====
Ben_Lat_plot =
  ggplot(data = filter(ABR_data_rat, current_inten == 80),
         aes(x = current_wave, y = Wave_lat, fill = genotype, group = interaction(current_wave, genotype, current_wave))) +
  geom_boxplot(linewidth = 1) +
  # geom_point(show.legend = FALSE, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "Wave",
       y = "Latancy @ 80dB",
       fill = "Genotype") +
  facet_wrap(~ current_freq, scales = "free_y") +
  theme_classic(base_size = 20) +
  theme(
    # legend.position=c(.8,.8)
    legend.position=c(.8,.2)
  )

print(Ben_Lat_plot)

ggsave(file="Fmr1_SD_ABR_lat.svg", path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben",
       plot = Ben_Lat_plot, 
       width = 10, height = 8, dpi = 600)

ggsave(file="Fmr1_SD_ABR_lat.jpg", path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben",
       plot = Ben_Lat_plot, 
       width = 10, height = 8, dpi = 600)

# Wave ratios -------------------------------------------------------------

W1_W4_ratio = 
ABR_data_rat %>%
  group_by(rat_ID, genotype, current_freq, , current_inten) %>%
  do(transmute(., current_wave = "1:4",
            Wave_amp = (filter(., current_wave == "1")$Wave_amp /
                        filter(., current_wave == "4")$Wave_amp),
            Wave_lat = (filter(., current_wave == "1")$Wave_lat /
                        filter(., current_wave == "4")$Wave_lat),
            RMS = filter(., current_wave == "1")$RMS) %>% unique)

W1_W5_ratio = 
  ABR_data_rat %>%
  group_by(rat_ID, genotype, current_freq, , current_inten) %>%
  do(transmute(., current_wave = "1:5",
               Wave_amp = (filter(., current_wave == "1")$Wave_amp /
                             filter(., current_wave == "5")$Wave_amp),
               Wave_lat = (filter(., current_wave == "1")$Wave_lat /
                             filter(., current_wave == "5")$Wave_lat),
               RMS = filter(., current_wave == "1")$RMS) %>% unique)

ABR_data_rat_complete = 
  mutate(ABR_data_rat, current_wave = as.character(current_wave)) %>%
  bind_rows(W1_W4_ratio, W1_W5_ratio) %>%
  mutate(current_wave = factor(current_wave, ordered = TRUE,
                               levels = c("1", "2", "3", "4", "5", "1:4", "1:5")))

## Graph ----
filter(ABR_data_rat_complete, current_freq == 0) %>%
  filter(current_wave %in% c("1", "1:5")) %>% 
  ggplot(aes(x = current_wave, y = Wave_lat, fill = genotype, group = interaction(current_wave, genotype, current_wave))) +
  geom_boxplot(linewidth = 1) +
  # geom_point(show.legend = FALSE, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "Wave",
       y = "BBN Amplitude",
       fill = "Genotype") +
  facet_wrap(~ current_inten, scales = "free_y") +
  theme_classic(base_size = 20) +
  theme(
    # legend.position=c(.8,.8)
    legend.position=c(.8,.1)
  )

# RMS ---------------------------------------------------------------------

## Graph ----
filter(ABR_data_rat_complete, current_freq == 0) %>%
  filter(current_wave %in% c("1")) %>% 
  ggplot(aes(x = current_inten, y = RMS, color = genotype, group = genotype)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 0, position = position_dodge(0.03)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(0.03), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.03), linewidth = 1)  +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_x_continuous(breaks = seq(from = 10, to = 90, by = 10)) +
  labs(x = "Intensity (dB)",
       y = "Signal Strength (RMS, \U003BCV)",
       color = "Genotype") +
  theme_classic(base_size = 20) +
  theme(
    # legend.position=c(.8,.8)
    legend.position=c(.8,.2)
  )
