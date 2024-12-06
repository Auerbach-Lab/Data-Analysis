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


# Piloting data checks ----------------------------------------------------

TH_table %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(Duration == 50) %>%
  filter(detail %in% c("Recheck", "Vehicle (Tween 80)", "Rapamycin (6mg/kg)", "None")) %>%
  reframe(TH = mean(TH, na.rm = TRUE),
          .by = c(genotype, detail))

vehicle_check_rxn =
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Recheck", "Vehicle (Tween 80)", "None")) %>%
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
  mutate(detail = str_replace(detail, pattern = "(Recheck|None)", replacement = "Baseline")) %>%
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
               geom = "line", position = position_dodge(0.5), linewidth = 1) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "Treatment") +
  scale_linetype_manual(values = c("Baseline" = "solid", "Vehicle (Tween 80)" = "dotted", 
                                   "Rapamycin (6mg/kg)" = "solid")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  theme_classic() +
  theme(
    legend.position = c(.8,.85),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 


# Rapamycin treatment -----------------------------------------------------
## Threshold -----

TH_table %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(Duration == 50) %>%
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)", "Recheck", "Post Treatment")) %>%
  reframe(TH = mean(TH, na.rm = TRUE),
          .by = c(genotype, detail)) %>%
  arrange(genotype, detail)

## Threshold Graph -----

TH_graph = 
TH_table %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(Duration == 50) %>%
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)",
                       "Recheck", "Post Treatment", "None")) %>%
  mutate(detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
         detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)",
                                            "Rapamycin (6mg/kg)", 
                                            "Post Treatment", "3+w Post Treatment",
                                            "Rapamycin 2 (6mg/kg)"),
                         labels = c("Baseline", "Vehicle", "Rapamycin",
                                    "Recovery", "3+ weeks Post Treatment", "Rapamycin 2"))) %>%
  ggplot(aes(x = genotype, y = TH, fill = genotype)) +
  geom_boxplot(linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  # stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
  #              position = position_dodge(1), vjust = 2, size = 3) +
  labs(x = "",
       y = "Threshold (dB)",
       fill = "Genotype") +
  scale_y_continuous(limits = c(19, 31), breaks = c(seq(18, 32, 2))) +
  scale_fill_manual(values = c("WT" = "darkgrey", "Het" = "deepskyblue", "KO" = "red")) +
  facet_wrap(~ detail, nrow = 1, strip.position = "top") +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(size = 18, face = "bold"),
    panel.grid.major.y = element_line(color = rgb(225, 225, 225, 255,
                                                  maxColorValue = 255))
  )

ggsave(filename = "TH_rapa.jpg",
       plot = TH_graph,
       width = 900, height = 900, units = "px", dpi = 100)

## Rxn data set -----
Rap_rxn =
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)", "Recheck",
                       "Post Treatment", "3+w Post Treatment", "Rapamycin 2 (6mg/kg)")) %>%
  filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  left_join(Tsc2_treatment_dates %>% 
              filter(detail == "Vehicle (Tween 80)") %>% 
              reframe(start_date = (date - 14) %>% str_remove_all(pattern = "-"), .by = rat_ID),
            by = join_by(rat_ID)) %>%
  rowwise() %>%
  filter(date > start_date) %>%
  unnest(reaction) %>% 
  reframe(Rxn = mean(Rxn, na.rm = TRUE) * 1000,
          .by = c(rat_ID, rat_name, sex, genotype, line, detail,
                  `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`)) %>%
  mutate(detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)",
                                            "Rapamycin (6mg/kg)", 
                                            "Post Treatment", "3+w Post Treatment",
                                            "Rapamycin 2 (6mg/kg)"),
                         labels = c("Baseline", "Vehicle", "Rapamycin",
                                    "Recovery", "3+ weeks Post Treatment", "Rapamycin 2")))

## Rxn Graph ====
Rap_rxn %>%
  filter(detail %in% c("Rapamycin", "Baseline")) %>%
  # filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  filter(`Inten (dB)` > 15 & `Inten (dB)` < 85) %>%
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
       color = "Genotype", linetype = "Treatment") +
  scale_linetype_manual(values = c("Baseline" = "solid", "Vehicle (Tween 80)" = "dotted", 
                                   "Rapamycin" = "longdash", "Post Treatment" = "dotdash",
                                   "3+w Post Treatment" = "twodash",
                                   "Rapamycin 2 (6mg/kg)" = 11)) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 


# WT graph ----------------------------------------------------------------
Rap_rxn %>%
  filter(genotype == "WT") %>%
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
       # title = glue("Wildtype"),
       color = "Genotype", linetype = "Treatment") +
  scale_linetype_manual(values = c("Recheck" = "longdash", "Vehicle (Tween 80)" = "dotted", 
                                   "Rapamycin (6mg/kg)" = "solid", "Post Treatment" = "dotdash",
                                   "3+w Post Treatment" = "twodash",
                                   "Rapamycin 2 (6mg/kg)" = 11)) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 


# Het graph ---------------------------------------------------------------
Rap_rxn %>%
  filter(genotype == "Het") %>%
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
       color = "Genotype", linetype = "Treatment") +
  scale_linetype_manual(values = c("Recheck" = "longdash", "Vehicle (Tween 80)" = "dotted", 
                                   "Rapamycin (6mg/kg)" = "solid", "Post Treatment" = "dotdash",
                                   "3+w Post Treatment" = "twodash",
                                   "Rapamycin 2 (6mg/kg)" = 11)) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 

# Individual Graphs -------------------------------------------------------

Individual_Graphs = 
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)", "Recheck",
                       "Post Treatment", "3+w Post Treatment",
                       "Rapamycin 2 (6mg/kg)", "Post Treatment 2", "None")) %>%
  filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  unnest(reaction) %>% 
  filter(`Inten (dB)` >= 15) %>%
  mutate(name = rat_name,
         detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
         detail = factor(detail, ordered = TRUE,
                         levels = c("Recheck", "Vehicle (Tween 80)", "Rapamycin (6mg/kg)",
                                    "Post Treatment", "3+w Post Treatment",
                                    "Rapamycin 2 (6mg/kg)", "Post Treatment 2"),
                         labels = c("Pre-Treatment", "Vehicle", "Rapamycin",
                                    "Recovery", "Post Treatment",
                                    "Rapamycin 2", "Post Treatment 2"))) %>%
  group_by(rat_ID, name) %>%
  do(single_rat_graph = 
       ggplot(data = .,
              aes(x = `Inten (dB)`, y = Rxn * 1000, color = detail)) +
       stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                    fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
                    fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
                    geom = "errorbar", position = position_dodge(0.5), width = 0) +
       stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                    geom = "point", position = position_dodge(0.5), size = 3) +
       stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                    geom = "line", position = position_dodge(0.5)) +
       # geom_smooth(se = FALSE, na.rm = TRUE) +
       labs(x = "Intensity (dB)",
            y = "Reaction time (ms, mean +/- SE)",
            color = "Treatment",
            title = glue("{unique(.$rat_name)} ({unique(.$sex)}, {unique(.$genotype)})")) +
       scale_color_manual(values = c("Pre-Treatment" = "black", 
                                     "Vehicle" = "darkblue", 
                                     "Rapamycin" = "red",
                                     "Recovery" = "goldenrod",
                                     "Post Treatment" = "forestgreen",
                                     "Rapamycin 2" = "deepskyblue",
                                     "Post Treatment 2" = "magenta")) +
       scale_x_continuous(breaks = seq(0, 90, by = 10)) +
       theme_classic() +
       theme(
         plot.title = element_text(hjust = 0.5),
         panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
         # legend.position = "none"
         # legend.position = c(.9,.85)
       )
  ) %>%
  arrange(name)

Individual_Graphs[c(1:8),]$single_rat_graph
