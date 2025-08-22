# Load Data (optional) ----------------------------------------------------

# source("Fmr1 vs Tsc2 in LE BBN tone data.R")

# Variables ---------------------------------------------------------------

save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Fmr1 LE"

# Graphing ----------------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}

n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}

# TH Averages -------------------------------------------------------------
FXS_TH_averages = TH_table %>%
  filter(line == "Fmr1-LE") %>%
  group_by(genotype, detail, Frequency , Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop") %>%
  arrange(Frequency, Duration, genotype)

print(FXS_TH_averages)

fwrite(FXS_TH_averages, file = glue("{save_folder}/Sample size.csv"))

# TH Graph ----------------------------------------------------------------

Fmr1_TH_gaph = 
  TH_table %>%
    filter(line == "Fmr1-LE") %>%
    filter(detail != "Rotating" & Frequency == 0) %>%
    {if (drop_TP3) filter(., rat_name != "TP3")} %>%
    filter(detail != "Rotating") %>%
    mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
    mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
             factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
    ggplot(aes(x = genotype, y = TH, shape = line,
               fill = genotype, color = group, group = interaction(group, line, genotype))) +
    geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
    # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
    stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
                 position = position_dodge(1), vjust = 2, size = 3) +
    # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
    scale_color_manual(values = c("Group 1" = "darkblue", "Group 2" = "goldenrod")) +
    scale_fill_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "lightcoral")) +
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

print(Fmr1_TH_gaph)

ggsave(filename = "Thresholds for BBN by group for Fmr1 LE.jpg",
       path = save_folder,
       plot = Fmr1_TH_gaph,
       width = 10, height = 8, units = "in", dpi = 300)

## Overall ----
TH_table %>%
  filter(line == "Fmr1-LE") %>%
  filter(Frequency == 0) %>%
  filter(detail != "Rotating") %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
           factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH,
             fill = genotype, group = genotype)) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
               position = position_dodge(1), vjust = 2, size = 3) +
  # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
  scale_fill_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "lightcoral")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       fill = "Genotype") +
  facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Thresholds for BBN for Fmr1 LE (summary fig).jpg",
       path = save_folder,
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 300)



# Reaction Graphs ------------------------------------------------------------

single_Frequency = "BBN"

Fmr1_single_frequency_rxn_graph = 
  Rxn_table %>%
    filter(line == "Fmr1-LE") %>%
    filter(Duration %in% c(50, 100, 300)) %>%
    # rename(Intensity = `Inten (dB)`) %>%
    filter(detail == "Alone") %>%
    mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
    # filter(rat_ID < 314) %>%
    mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
    filter(! str_detect(Intensity, pattern = "5$")) %>%
    filter(Intensity < 90 & Intensity > 10) %>%
    filter(Frequency == single_Frequency) %>%
    ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(group),
               color = genotype, group = interaction(Duration, group, genotype))) +
    ## Overall average lines
    stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
                 fun = function(x) mean(x, na.rm = TRUE),
                 fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
                 fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
                 geom = "errorbar", width = 1.5, linetype = "solid", linewidth = 1.5, alpha = 0.5,
                 position = position_dodge(4)) +
    stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
                 fun = function(x) mean(x, na.rm = TRUE),
                 geom = "line", linetype = "solid", linewidth = 1.5, alpha = 0.5,
                 position = position_dodge(4)) +
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
    scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
    scale_x_continuous(breaks = seq(0, 90, by = 10)) +
    facet_wrap(~ Duration) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    ) 

print(Fmr1_single_frequency_rxn_graph)

ggsave(filename = "Reaction times BBN durations by group.jpg",
       path = save_folder,
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 300)
## Temportal Integration ----
Rxn_table %>%
  filter(line == "Fmr1-LE") %>%
  filter(Duration %in% c(50, 100, 300)) %>%
  # rename(Intensity = `Inten (dB)`) %>%
  filter(detail == "Alone") %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  # filter(rat_ID < 314) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == single_Frequency) %>%
  ggplot(aes(x = Intensity, y = Rxn, 
             # linetype = as.factor(Duration), shape = as.factor(Duration),
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
       color = "Genotype", linetype = "Duration", shape = "Duration") +
  # scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

ggsave(filename = "Reaction times BBN durations (summary figure).jpg",
       path = save_folder,
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 300)

## Temportal Integration ----
Rxn_table %>%
  filter(line == "Fmr1-LE") %>%
  filter(Duration %in% c(50, 100, 300)) %>%
  # rename(Intensity = `Inten (dB)`) %>%
  filter(detail == "Alone") %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  # filter(rat_ID < 314) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == single_Frequency) %>%
  ggplot(aes(x = Intensity, y = Rxn, 
             linetype = as.factor(Duration), shape = as.factor(Duration),
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
       color = "Genotype", linetype = "Duration", shape = "Duration") +
  # scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ genotype) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

ggsave(filename = "Reaction times BBN durations (temporal integration figure).jpg",
       path = save_folder,
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 300)


# Individual Graphs -------------------------------------------------------

BBN_Individual_Graphs = 
  Rxn_table %>%
  filter(line == "Fmr1-LE") %>%
  # filter(Duration %in% c(300, 100)) %>%
  filter(detail == "Alone") %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  # filter(rat_ID < 314) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == "BBN") %>%
  mutate(name = rat_name) %>%
  group_by(rat_ID, name) %>%
  do(bbn_single_rat_graph = 
       ggplot(data = .,
              aes(x = Intensity, y = Rxn, linetype = as.factor(Duration),
                  color = genotype, group = interaction(Duration, group, genotype))) +
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
            color = "Genotype", linetype = "",
            title = glue("{unique(.$rat_name)}")) +
       scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
       scale_x_continuous(breaks = seq(0, 90, by = 10)) +
       theme_classic() +
       theme(
         plot.title = element_text(hjust = 0.5),
         panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
       )
  ) %>%
  arrange(rat_ID)

BBN_Individual_Graphs$bbn_single_rat_graph
