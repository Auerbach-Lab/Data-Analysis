ggplot(discrimination_FA_table_by_type %>%
         filter(detail == "Normal") %>%
         filter(line == "Tsc2") %>% filter(Type != "Broad"), 
       aes(x = octave_steps, y = FA_percent_detailed * 100,
           color = genotype, fill = Type, 
           group = interaction(Type, detail, line, genotype))) +
  geom_hline(yintercept = 50, color = "forestgreen", linewidth = 1.5) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.5, linewidth = 1.0, position = position_dodge(0.2)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
  # mean for each frequency by genotype
  stat_summary(aes(shape = Type), fun = mean, geom = "point",
               position = position_dodge(.2), size = 2, stroke = 2) +
  # add labels by x-axis
  geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
                          tone = "No Go", genotype = "WT", line = "Tsc2", detail = "Normal", Type = "Broad"),
            aes(label = tone), size = 4, show.legend = FALSE, vjust = 2, hjust = -0.2,
            family = "EconSansCndReg") +
  # geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
  #                         tone = "Go", genotype = "WT", line = "Tsc2", detail = "Normal", Type = "Zoom"),
  #           aes(label = tone), size = 4, show.legend = FALSE, vjust = 3.1, hjust = 4,
  #           family = "EconSansCndReg") +
  # table on graph
  annotate(geom = "table", x = 12, y = 100,
           label = list(discrimination_FA_n_table %>%
                          filter(detail %in% c("Normal")) %>%
                          filter(str_detect(Genotype, pattern = "Tsc2")) %>%
                          select(Genotype, n)),
           table.theme = ttheme_gtplain(
             padding = unit(c(1, 0.75), "char")
           )) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_shape_manual(values = c("Broad" = 21, "Zoom" = 24)) +
  scale_fill_manual(values = c("Broad" = "slategrey", "Zoom" = "tan4")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       title = "Discrimination across an octave",
       fill = "Octave Steps", shape = "Octave Steps",
       color = "Genotype", linetype = "Octave Steps") +
  theme_ipsum_es()

ggsave(filename = "Tsc2 octave discrimination - Zoom.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 11, height = 8, units = "in", dpi = 150)


ggplot(discrimination_FA_table %>%
         filter(detail == "Normal") %>%
         filter(line == "Tsc2"), 
       aes(x = octave_steps, y = FA_percent_detailed * 100,
           color = genotype, 
           group = interaction(detail, line, genotype))) +
  geom_hline(yintercept = 50, color = "forestgreen", linewidth = 1.5) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0.5, linewidth = 1.0, position = position_dodge(0.2)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
  # mean for each frequency by genotype
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(.2), size = 2, stroke = 2) +
  # add labels by x-axis
  geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
                          tone = "No Go", genotype = "WT", line = "Tsc2", detail = "Normal", Type = "Broad"),
            aes(label = tone), size = 4, show.legend = FALSE, vjust = 2, hjust = -0.2,
            family = "EconSansCndReg") +
  # geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
  #                         tone = "Go", genotype = "WT", line = "Tsc2", detail = "Normal", Type = "Zoom"),
  #           aes(label = tone), size = 4, show.legend = FALSE, vjust = 3.1, hjust = 4,
  #           family = "EconSansCndReg") +
  # table on graph
  annotate(geom = "table", x = 12, y = 100,
           label = list(discrimination_FA_n_table %>%
                          filter(detail %in% c("Normal")) %>%
                          filter(str_detect(Genotype, pattern = "Tsc2")) %>%
                          select(Genotype, n)),
           table.theme = ttheme_gtplain(
             padding = unit(c(1, 0.75), "char")
           )) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_shape_manual(values = c("Broad" = 21, "Zoom" = 24)) +
  scale_fill_manual(values = c("Broad" = "slategrey", "Zoom" = "tan4")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "blue", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       title = "Discrimination across an octave",
       fill = "Octave Steps", shape = "Octave Steps",
       color = "Genotype", linetype = "Octave Steps") +
  theme_ipsum_es()

ggsave(filename = "Tsc2 octave discrimination - combined.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 11, height = 8, units = "in", dpi = 150)
