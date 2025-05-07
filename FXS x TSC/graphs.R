# Load Data (optional) ----------------------------------------------------

## Load data ---
source("data.R")

## Process data ---
# source("data.R")

# Graphing ----------------------------------------------------------------
n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}


# Get Group info ----------------------------------------------------------

## Group 1 ----
Group1 = # c(745, 746, 747, 748)
  # dynamically selected based on start date from October 2024
  core_data %>%
  filter(date > 20241001 & date < 20241201) %>% 
  .$rat_ID %>% # use rat_ID because its unique
  unique # de-duplicate

## Group 2 ----
Group2 = # c(769, 770, 772, 774) & c(744, 749, 771, 775)
  # dynamically selected based on start date of January 2025
  core_data %>%
  filter(date > 20250101 & date < 20250112) %>% 
  filter(! rat_ID %in% Group1) %>%
  .$rat_ID %>% # use rat_ID because its unique
  unique # de-duplicate

# Individual Graphs -------------------------------------------------------

Individual_Graphs = 
  core_data %>%
  filter(! task %in% c("Reset")) %>%    # Omit Training & Reset days
  filter(! task %in% c("Training")) %>%    # Omit Training & Reset days
  filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  unnest(reaction) %>% 
    filter(! (stim_type == "tone" & `Dur (ms)` %in% c(300, 100))) %>%    # Only have 50ms for tone data
  mutate(name = rat_name) %>%
  group_by(rat_ID, name) %>%
  do(single_rat_graph = 
       ggplot(data = .,
              aes(x = `Inten (dB)`, y = Rxn * 1000, 
                  color = as.factor(`Dur (ms)`), linetype = detail, 
                  shape = as.factor(`Freq (kHz)`),
                  group = interaction(as.factor(`Dur (ms)`), detail, as.factor(`Freq (kHz)`)))) +
       stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                    fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
                    fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
                    geom = "errorbar", position = position_dodge(0.5), width = 0) +
       stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                    geom = "point", position = position_dodge(0.5), size = 2) +
       stat_summary(fun = function(x) mean(x, na.rm = TRUE),
                    geom = "line", position = position_dodge(0.5)) +
       scale_linetype_manual(values = c("None" = "solid", "Recheck" = "dashed", "Mixed" = "solid")) +
       # geom_smooth(se = FALSE, na.rm = TRUE) +
       labs(x = "Intensity (dB)",
            y = "Reaction time (ms, mean +/- SE)",
            color = "Stim Duration",
            linetype = "Exp. Detail",
            shape = "Frequency",
            title = glue("{unique(.$rat_name)} ({unique(.$genotype)})")) +
       facet_wrap(~ stim_type, scales = "free_y", nrow = 2) + 
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

Individual_Graphs %>%
  # filter(rat_ID %in% Group1) %>%
  filter(name %in% c("Orange2", "Orange3", "Orange4", #Orange1, #(Group 1.1), 
                     "Lime1", "Lime2", "Orange5", "Orange6",  #still in training (Group 1.2)
                     "Blue1", "Blue2", "Purple6", "Lime6")) %>% #(Group 1.3)
  .$single_rat_graph

# TH Graph ----------------------------------------------------------------

## BBN TH by Genotype ----
TH_table %>%
  filter(Frequency == 0) %>%
  mutate(Duration = as.factor(Duration)) %>%
  ggplot(aes(x = Duration, y = TH,
             fill = genotype, group = interaction(Duration, genotype))) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
               position = position_dodge(1), vjust = 2, size = 3) +
  scale_fill_manual(values = c("Wild-type" = "black", "Double KO" = "darkmagenta",
                               "TSC only" = "deepskyblue", "FXS only" = "red")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       fill = "Genotype") +
  facet_wrap( ~ genotype, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## BBN TH by Duration ----
TH_table %>%
  mutate(Duration = as.factor(Duration)) %>%
  ggplot(aes(x = genotype, y = TH,
             fill = genotype, group = genotype)) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
               position = position_dodge(1), vjust = 2, size = 3) +
  scale_fill_manual(values = c("Wild-type" = "black", "Double KO" = "darkmagenta",
                                "TSC only" = "deepskyblue", "FXS only" = "red")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       fill = "Genotype") +
  facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


## Tones TH by Genotype ----
TH_table %>%
  filter(Duration == 50) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  ggplot(aes(x = Frequency, y = TH,
             fill = genotype, group = interaction(Frequency, genotype))) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
               position = position_dodge(1), vjust = 2, size = 3) +
  scale_fill_manual(values = c("Wild-type" = "black", "Double KO" = "darkmagenta",
                               "TSC only" = "deepskyblue", "FXS only" = "red")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       fill = "Genotype") +
  facet_wrap( ~ genotype, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

# Reaction time -----------------------------------------------------------

## BBN ----
Rxn_table %>%
  filter(Frequency == 0) %>%
  filter(Intensity >= 20) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  # filter(! str_detect(Intensity, pattern = "5$")) %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(Duration),
             color = genotype, group = interaction(Duration, genotype))) +
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
       color = "Genotype", linetype = "") +
  scale_color_manual(values = c("Wild-type" = "black", "Double KO" = "darkmagenta",
                                "TSC only" = "deepskyblue", "FXS only" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

## Temporal Integration ----
Rxn_table %>%
  filter(Frequency == 0) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(Intensity >= 40) %>%
  # filter(! str_detect(Intensity, pattern = "5$")) %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(Duration),
             color = genotype, group = interaction(Duration, genotype))) +
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
       color = "Genotype", linetype = "") +
  scale_color_manual(values = c("Wild-type" = "black", "Double KO" = "darkmagenta",
                                "TSC only" = "deepskyblue", "FXS only" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ genotype) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 


## Tones ----
Rxn_table %>%
  filter(Duration == 50) %>%
  filter(Intensity >= 20) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  # filter(! str_detect(Intensity, pattern = "5$")) %>%
  ggplot(aes(x = Intensity, y = Rxn,
             color = genotype, group = interaction(Frequency, genotype))) +
  ## Lines for each group
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 1.5, position = position_dodge(0.5)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(0.5), size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line", position = position_dodge(0.5)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "") +
  scale_color_manual(values = c("Wild-type" = "black", "Double KO" = "darkmagenta",
                                "TSC only" = "deepskyblue", "FXS only" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Frequency, ncol = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 
