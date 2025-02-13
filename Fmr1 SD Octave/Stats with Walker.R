
# Data --------------------------------------------------------------------
walker_data = fread(glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Walker/Fmr1_SD_Octave_20240426.csv"))
walker_data1 = fread(glue("C:/Users/Noelle/Box/tone_discrim_2024/final_docs/csvs_for_stats_check/octave_stats_plot_works.csv"))

# FA% ---------------------------------------------------------------------
## Data ----
FA_stats_data =
Discrimination_data %>%
  filter(detail == "Normal") %>%
  # filter(Type == "Zoom") %>%
  # filter(octave_steps != 9) %>%
  group_by(rat_ID, rat_name, Genotype, octave_steps) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")


FA_stats_data_walker =
walker_data1 %>%
  filter(detail == "Normal") %>%
  filter(HL_state == "baseline") %>%
  # filter(Type == "Zoom") %>%
  # filter(octave_steps != 9) %>%
  group_by(rat_ID, rat_name, Genotype, octave_steps) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")

## Graph ----
ggplot(data = FA_stats_data_walker,
       aes(x = octave_steps, y = FA_percent_detailed * 100,
           color = Genotype, group = interaction(Genotype))) +
  geom_hline(yintercept = 50, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  # geom_smooth(se = FALSE) +
  # add the TH line
  # geom_vline(data = TH %>% filter(detail == "Normal"), aes(xintercept = TH, color = Genotype), show.legend = FALSE) +
  # mean for each frequency by genotype
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  # add labels by x-axis
  # geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
  #                         tone = "No Go", Genotype = "WT"),
  #           aes(label = tone), size = 5, show.legend = FALSE, vjust = 1,
  #           family = "EconSansCndReg") +
  # geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
  #                         tone = "Go", Genotype = "WT"),
  #           aes(label = tone), size = 4, show.legend = FALSE, vjust = 1,
  #           family = "EconSansCndReg") +
  scale_x_continuous(breaks = c(1, seq(0, 12, by = 2))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       # title = "Discrimination across an octave",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

### Save ----
# ggsave(file="FA_combined.svg",
#        path = glue(Box_folder, 
#                    "/tone_discrim_2024/final_docs"),
#        plot = last_plot(),
#        width=10, height=8)


## Stats ----

FA_aov = 
  aov(FA_percent_detailed ~ Genotype * octave_steps,
      data = FA_stats_data_walker %>%
        filter(octave_steps != 9)  %>%
        mutate(octave_steps = ordered(octave_steps)))

### Normality ----
shapiro.test(FA_aov$residuals)$p.value

summary(FA_aov)

broom::tidy(TukeyHSD(FA_aov)) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  # filter(sig != " ") %>%
  mutate(Comp1 = str_split_fixed(.$contrast, '-', 2)[,1],
         Comp2 = str_split_fixed(.$contrast, '-', 2)[,2],
         geno1 = str_split_fixed(Comp1, ':', 3)[,1],
         octave_step1 = str_split_fixed(Comp1, ':', 3)[,2] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, ':', 3)[,1],
         octave_step2 = str_split_fixed(Comp2, ':', 3)[,2] %>% as.numeric(),) %>%
  filter(geno1 != geno2) %>%
  filter(octave_step1 == octave_step2) %>%
  View



# Broad vs Zoom main effect -----------------------------------------------

aov(FA_percent_detailed ~ Genotype * octave_steps * Type,
    data =   Discrimination_data %>%
      filter(detail == "Normal") %>%
      filter(octave_steps != 9) %>%
      group_by(rat_ID, rat_name, Genotype, octave_steps, Type) %>%
      summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
                .groups = "drop") %>%
      filter(octave_steps != 9)  %>%
      mutate(octave_steps = ordered(octave_steps))) %>%
  summary()

aov(FA_percent_detailed ~ Genotype * octave_steps,
    data =   Discrimination_data %>%
      filter(detail == "Normal") %>%
      filter(Type == "Zoom") %>%
      filter(octave_steps %in% c(5, 6, 7)) %>%
      group_by(rat_ID, rat_name, Genotype, octave_steps, Type) %>%
      summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
                .groups = "drop") %>%
      filter(octave_steps != 9)  %>%
      mutate(octave_steps = ordered(octave_steps))) %>%
  summary()

aov(FA_percent_detailed ~ Genotype * octave_steps,
    data =   Discrimination_data %>%
      filter(detail == "Normal") %>%
      filter(Type == "Broad") %>%
      # filter(octave_steps != 9) %>%
      group_by(rat_ID, rat_name, Genotype, octave_steps, Type) %>%
      summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
                .groups = "drop") %>%
      filter(octave_steps != 9)  %>%
      mutate(octave_steps = ordered(octave_steps))) %>%
  summary()

# FA by Zoom --------------------------------------------------------------

## Data ----
FA_zoom_stats_data =
  Discrimination_data %>%
  filter(detail == "Normal") %>%
  # filter(Type == "Zoom") %>%
  # filter(octave_steps != 9) %>%
  group_by(rat_ID, rat_name, Genotype, octave_steps, Type) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")

## graph ----
ggplot(data = FA_zoom_stats_data %>% filter(Type == "Zoom"),
       aes(x = octave_steps, y = FA_percent_detailed * 100,
           color = Genotype, group = interaction(Genotype))) +
  geom_hline(yintercept = 50, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  # geom_smooth(se = FALSE) +
  # add the TH line
  # geom_vline(data = TH %>% filter(detail == "Normal"), aes(xintercept = TH, color = Genotype), show.legend = FALSE) +
  # mean for each frequency by genotype
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  # add labels by x-axis
  # geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
  #                         tone = "No Go", Genotype = "WT"),
  #           aes(label = tone), size = 5, show.legend = FALSE, vjust = 1,
  #           family = "EconSansCndReg") +
  # geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
  #                         tone = "Go", Genotype = "WT"),
  #           aes(label = tone), size = 4, show.legend = FALSE, vjust = 1,
  #           family = "EconSansCndReg") +
  scale_x_continuous(breaks = c(1, seq(0, 12, by = 2))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       # title = "Discrimination across an octave",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))


## Stats ----

FA_Zoom_aov = 
  aov(FA_percent_detailed ~ Genotype * octave_steps * Type,
      data = FA_zoom_stats_data %>%
        filter(octave_steps != 9)  %>%
        mutate(octave_steps = ordered(octave_steps)))

### Normality ----
shapiro.test(FA_Zoom_aov$residuals)$p.value

summary(FA_Zoom_aov)

broom::tidy(TukeyHSD(FA_Zoom_aov)) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  # filter(sig != " ") %>%
  mutate(Comp1 = str_split_fixed(.$contrast, '-', 2)[,1],
         Comp2 = str_split_fixed(.$contrast, '-', 2)[,2],
         geno1 = str_split_fixed(Comp1, ':', 3)[,1],
         octave_step1 = str_split_fixed(Comp1, ':', 3)[,2] %>% as.numeric(),
         type1 = str_split_fixed(Comp1, ':', 3)[,3],
         geno2 = str_split_fixed(Comp2, ':', 3)[,1],
         octave_step2 = str_split_fixed(Comp2, ':', 3)[,2] %>% as.numeric(),
         type2 = str_split_fixed(Comp2, ':', 3)[,3],) %>%
  filter(geno1 != geno2) %>%
  filter(octave_step1 == octave_step2) %>%
  filter(! is.na(estimate)) %>%
  View

### step 5 only ----
FA_Zoom_5_aov = 
  aov(FA_percent_detailed ~ Genotype,
      data = FA_zoom_stats_data %>%
        filter(octave_steps == 5))

shapiro.test(FA_Zoom_5_aov$residuals)$p.value

summary(FA_Zoom_5_aov)

t.test(FA_percent_detailed ~ Genotype,
    data = FA_zoom_stats_data %>%
      filter(octave_steps == 7))



# Binned analysis ---------------------------------------------------------
FA_bin_stats_data =
  walker_data1 %>%
  filter(detail == "Normal") %>%
  filter(Type == "Zoom") %>%
  # filter(octave_steps != 9) %>%
  mutate(bin = case_when(octave_steps %in% c(1, 2, 3, 4) ~ "Low",
                         octave_steps %in% c(5, 6, 7, 8) ~ "Mid",
                         octave_steps %in% c(9, 10, 11, 12) ~ "High") %>%
           factor( levels = c("Low", "Mid", "High"))) %>%
  # filter(bin %in% c("Low", "Mid")) %>%
  group_by(rat_ID, rat_name, Genotype, bin) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")

## Graph ----
FA_bin_stats_data %>%
  ggplot(aes(x = bin, y = FA_percent_detailed, fill = Genotype, color = Genotype,
             group = interaction(bin, Genotype))) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ Type) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

### Save ----
# ggsave(file="FA_binned.svg",
#        path = glue(Box_folder,
#                    "/tone_discrim_2024/final_docs"),
#        plot = last_plot(),
#        width=10, height=8)

## Stats ----
FA_bin_aov =
aov(FA_percent_detailed ~ Genotype * bin,
    data = FA_bin_stats_data)

shapiro.test(FA_bin_aov$residuals)$p.value

summary(FA_bin_aov)

broom::tidy(TukeyHSD(FA_bin_aov)) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  # filter(sig != " ") %>%
  mutate(Comp1 = str_split_fixed(.$contrast, '-', 2)[,1],
         Comp2 = str_split_fixed(.$contrast, '-', 2)[,2],
         geno1 = str_split_fixed(Comp1, ':', 3)[,1],
         octave_step1 = str_split_fixed(Comp1, ':', 3)[,2],
         geno2 = str_split_fixed(Comp2, ':', 3)[,1],
         octave_step2 = str_split_fixed(Comp2, ':', 3)[,2]) %>%
  filter(geno1 != geno2) %>%
  filter(octave_step1 == octave_step2)

t.test(FA_percent_detailed ~ Genotype,
       data = FA_bin_stats_data %>%
         filter(bin == "Mid"))


# ABR thresholds ----------------------------------------------------------

walker_ABR_data = fread(glue(Box_folder, 
                             "/tone_discrim_2024/final_docs/csvs_for_stats_check/abr_thresholds.csv")) %>%
  as_tibble() 


## Graph ---
walker_ABR_data %>%
  summarise(TH = mean(current_inten, na.rm = TRUE),
            .by = c(rat_ID, genotype, current_freq)) %>%
  ggplot(aes(x = genotype, y = TH, fill = genotype, group = interaction(genotype))) +
  geom_boxplot() +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  facet_wrap(~ current_freq) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Stats ----

t.test(TH ~ genotype,
       data =  walker_ABR_data %>%
         summarise(TH = mean(current_inten, na.rm = TRUE),
                   .by = c(rat_ID, genotype)))


TH_aov =
  aov(TH ~ genotype * current_freq,
      data = walker_ABR_data %>%
        summarise(TH = mean(current_inten, na.rm = TRUE),
                  .by = c(rat_ID, genotype, current_freq)))

shapiro.test(TH_aov$residuals)$p.value

summary(TH_aov)


# d' TH -------------------------------------------------------------------

walker_dprime_table = fread(glue(Box_folder, 
                             "/tone_discrim_2024/final_docs/csvs_for_stats_check/",
                             "dprime_threshold_detection.csv")) %>%
  as_tibble() 

walker_dprime_data = 
  walker_dprime_table %>%
  mutate(Freq = str_extract(file_name, pattern = "(([:digit:]|-)+?kHz|BBN)"),
         genotype = str_extract(Genotype, pattern = "(WT|KO)$")) %>%
  summarise(dprime = mean(dprime, na.rm = TRUE),
            .by = c(rat_name, genotype, Freq, dB))

## graph ----

walker_dprime_data %>%
ggplot(aes(x = dB, y = dprime,
           color = genotype, group = interaction(genotype))) +
  geom_hline(yintercept = 1.5, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  scale_x_continuous(breaks = c(seq(-50, 90, by = 10))) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  facet_wrap(~ Freq) +
  labs(x = "Intensity",
       y = "d'",
       # title = "Discrimination across an octave",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))


## Stats ----
dprime_aov =
  aov(dprime ~ genotype * Freq,
      data = walker_dprime_data)

shapiro.test(dprime_aov$residuals)$p.value

summary(TH_aov)

# Kruskal Testing - Main effects only 
kruskal.test(dprime ~ genotype,
             data = walker_dprime_data %>%
               filter(Freq == "4-32kHz")) %>%
  tidy()


# Hits by day -------------------------------------------------------
walker_by_day_data = fread(glue(Box_folder, 
                                 "/tone_discrim_2024/final_docs/csvs_for_stats_check/updated_versions/",
                                 "fmr1_octave_over_time.csv")) %>%
  as_tibble() %>%
  mutate(day = as.numeric(day))

## Graph ----
walker_by_day_data %>%
  filter(type == "HIT") %>%
  ggplot(aes(x = day, y = rate, color = genotype)) +
  stat_summary(aes(group = rat),
               fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1, position = position_dodge(.1), alpha = 0.2) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  ylim(60, 100) +
  labs(x = "Day of training",
       y = "Hit %",
       # title = "Discrimination across an octave",
       color = "Genotype") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Stats ----

wilcox.test(hit ~ genotype,
       data =  walker_by_day_data %>%
         filter(type == "HIT") %>%
         summarise(hit = mean(rate, na.rm = TRUE),
                   .by = c(rat, genotype)))


Hit_by_Freq_aov =
  aov(Gaus ~ genotype * day,
      data = walker_by_day_data %>%
        filter(type == "HIT") %>%
        mutate(Gaus = LambertW::Gaussianize(rate)))

shapiro.test(Hit_by_Freq_aov$residuals)$p.value

summary(Hit_by_Freq_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("genotype", "day", # Main effects
         "interaction(genotype, day)"
), 
function(x) kruskal.test(reformulate(x, "rate"),
                         data = walker_by_day_data %>%
                           filter(type == "HIT"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)


# FA by day ---------------------------------------------------------------

## Graph ----
walker_by_day_data %>%
  filter(type == "FA") %>%
  ggplot(aes(x = day, y = rate, color = genotype)) +
  stat_summary(aes(group = rat),
               fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1, position = position_dodge(.1), alpha = 0.2) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  ylim(0, 100) +
  labs(x = "Day of training",
       y = "False Alarm %",
       # title = "Discrimination across an octave",
       color = "Genotype") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Stats ----

wilcox.test(hit ~ genotype,
       data =  walker_by_day_data %>%
         filter(type == "FA") %>%
         summarise(hit = mean(rate, na.rm = TRUE),
                   .by = c(rat, genotype)))


FA_by_date_aov =
  aov(Gaus ~ genotype * day,
      data = walker_by_day_data %>%
        filter(type == "FA") %>%
        mutate(Gaus = LambertW::Gaussianize(rate)))

shapiro.test(FA_by_date_aov$residuals)$p.value

summary(FA_by_date_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("genotype", "day" # Main effects
), 
function(x) kruskal.test(reformulate(x, "rate"),
                         data = walker_by_day_data %>%
                           filter(type == "FA"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)


### Post-Hoc Dunn's Test ----
FSA::dunnTest(rate ~ interaction(genotype, day),
              data = walker_by_day_data %>%
                filter(type == "FA"),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         date1 = str_split_fixed(Comp1, '\\.', 2)[,2] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         date2 = str_split_fixed(Comp2, '\\.', 2)[,2] %>% as.numeric()) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(date1 == date2) %>%
  arrange(date1) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2, -date1)



# Behavior TH -------------------------------------------------------------

walker_behavior_TH_table = fread(glue(Box_folder, 
                                      "/tone_discrim_2024/final_docs/csvs_for_stats_check/",
                                      "fmr1_sd_thresholds_bda.csv")) %>%
  as_tibble() 


## Graph ---
walker_behavior_TH_table %>%
  mutate(freq = as.factor(freq)) %>%
  ggplot(aes(x = freq, y = thresh, fill = genotype, group = interaction(genotype, freq))) +
  geom_boxplot() +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  # scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ freq) +
  labs(x = "Genotype",
       y = "Threshold",
       # title = "Discrimination across an octave",
       color = "Genotype") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Stats ----

wilcox.test(thresh ~ genotype,
       data =  walker_behavior_TH_table)


TH_behavior_aov =
  aov(thresh ~ genotype * freq,
      data = walker_behavior_TH_table%>%
        mutate(Gaus = LambertW::Gaussianize(thresh)))

shapiro.test(TH_behavior_aov$residuals)$p.value

summary(TH_behavior_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("genotype", "freq" # Main effects
), 
function(x) kruskal.test(reformulate(x, "thresh"),
                         data = walker_behavior_TH_table)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

### Post-Hoc Dunn's Test ----
FSA::dunnTest(thresh ~ interaction(genotype, freq),
              data = walker_behavior_TH_table,
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         date1 = str_split_fixed(Comp1, '\\.', 2)[,2] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         date2 = str_split_fixed(Comp2, '\\.', 2)[,2] %>% as.numeric()) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(date1 == date2) %>%
  arrange(date1) %>%
  rename(Frequency = date2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2, -date1)


# AC modeling ----------------------------------------------------------------

acx_model = fread(glue(code_folder,
                       "acx_tone_io_params_stats.csv")) %>%
  as_tibble()


## Graph by channel ----
acx_model %>%
  ggplot(aes(x = channel, y = slope, 
             fill = genotype, color = genotype, group = interaction(channel, genotype))) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ Type) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Graph by file ----
acx_model %>%
  ggplot(aes(x = file, y = slope, 
             fill = genotype, color = genotype, group = interaction(file, genotype))) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ Type) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## mixed graph ----
acx_model %>%
  pivot_longer(c(slope:y_max)) %>%
  ggplot(aes(x = genotype, y = value, 
             fill = genotype, color = genotype)) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  facet_wrap(~ name, scale = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


## Stats ----
wilcox.test(slope ~ genotype,
            data =  acx_model)


acx_model_aov =
  aov(slope ~ genotype * channel,
      data = acx_model %>%
        mutate(Gaus = LambertW::Gaussianize(slope)))

shapiro.test(acx_model_aov$residuals)$p.value

summary(acx_model_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("slope", "threshold", "y_min", "y_max" # Main effects
), 
function(x) kruskal.test(reformulate("genotype", x),
                         data = acx_model)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

### Post-Hoc Dunn's Test ----
FSA::dunnTest(slope ~ interaction(genotype, file),
              data = acx_model,
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         date1 = str_split_fixed(Comp1, '\\.', 2)[,2] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         date2 = str_split_fixed(Comp2, '\\.', 2)[,2] %>% as.numeric()) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(date1 == date2) %>%
  arrange(date1) %>%
  rename(Channel = date2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2, -date1)

# AC Gaus modeling ----------------------------------------------------------------
acx__gaus_model = fread(glue(Box_folder,
                       "/tone_discrim_2024/",
                       "acx_gauss_paramsfra.csv")) %>%
  as_tibble() %>%
  rename(genotype = Genotype)


## Graph ----
acx__gaus_model %>%
  ggplot(aes(x = channel, y = gain, 
             fill = genotype, color = genotype, group = interaction(channel, genotype))) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ Type) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## mixed graph ----
acx__gaus_model %>%
  filter(y_max < 150) %>%
  pivot_longer(c(y_min:y_max, gain, range)) %>%
  ggplot(aes(x = genotype, y = value, 
             fill = genotype, color = genotype)) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  facet_wrap(~ name, scale = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


## Stats ----
wilcox.test(gain ~ genotype,
            data =  acx__gaus_model)


acx_model_gaus_aov =
  aov(gain ~ genotype * channel,
      data = acx__gaus_model %>%
        mutate(Gaus = LambertW::Gaussianize(gain)))

shapiro.test(acx_model_gaus_aov$residuals)$p.value

summary(acx_model_gaus_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("gain", "range", "y_min", "y_max" # Main effects
), 
function(x) kruskal.test(reformulate("genotype", x),
                         data = acx__gaus_model)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

### Post-Hoc Dunn's Test ----
FSA::dunnTest(slope ~ interaction(genotype, channel),
              data = acx__gaus_model,
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
         date1 = str_split_fixed(Comp1, '\\.', 2)[,2] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
         date2 = str_split_fixed(Comp2, '\\.', 2)[,2] %>% as.numeric()) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(date1 == date2) %>%
  arrange(date1) %>%
  rename(Channel = date2) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2, -date1)

# IC Gaus modeling ----------------------------------------------------------------

ic_model_gaus = fread(glue(Box_folder,
                       "/tone_discrim_2024/",
                       "ic_gauss_paramsfra.csv")) %>%
  as_tibble() %>%
  rename(genotype = Genotype)


## Graph ----
ic_model_gaus %>%
  ggplot(aes(x = channel, y = gain, 
             fill = genotype, color = genotype, group = interaction(channel, genotype))) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ Type) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## mixed graph ----
ic_model_gaus %>%
  pivot_longer(c(y_min:y_max, gain, range)) %>%
  ggplot(aes(x = genotype, y = value, 
             fill = genotype, color = genotype)) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  facet_wrap(~ name, scale = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


## Stats ----
ic_model_gaus_aov =
  aov(gain ~ genotype * channel,
      data = ic_model_gaus %>%
        mutate(Gaus = LambertW::Gaussianize(gain)))

shapiro.test(ic_model_gaus_aov$residuals)$p.value

summary(ic_model_gaus_aov)


# IC modeling ----------------------------------------------------------------

ic_model = fread(glue(Box_folder,
                      "/tone_discrim_2024/",
                      "ic_params_io_rlf_final.csv")) %>%
  as_tibble()


## Graph ----
ic_model %>%
  ggplot(aes(x = channel, y = slope, fill = genotype, color = genotype, 
             group = interaction(channel, genotype))) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  # facet_wrap(~ Type) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## mixed graph ----
ic_model %>%
  pivot_longer(c(slope:y_max)) %>%
  ggplot(aes(x = genotype, y = value, 
             fill = genotype, color = genotype)) +
  geom_boxplot(color = "black") +
  geom_point(position = position_dodge(0.75)) +
  scale_color_manual(values = c("WT" = "blue", "KO" = "violet")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  facet_wrap(~ name, scale = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


## Stats ----
wilcox.test(slope ~ genotype,
            data =  ic_model)


IC_model_aov =
  aov(slope ~ genotype * channel,
      data = ic_model %>%
        mutate(Gaus = LambertW::Gaussianize(slope)))

shapiro.test(IC_model_aov$residuals)$p.value

summary(IC_model_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("slope", "threshold", "y_min", "y_max" # Main effects
), 
function(x) kruskal.test(reformulate("genotype", x),
                         data = acx_model)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)


