# Graphing ----------------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}



# Task grouping -----------------------------------------------------------

Rxn_table %>%
  filter(task %in% c("Rotating", "CNO 3mg/kg", "Saline")) %>%
  filter(rat_ID != 100) %>%
  ggplot(aes(x = Position, y = Rxn, color = task)) +
  # geom_point(aes(group = ID, color = Genotype), alpha = 0.3)+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 1, position = position_dodge(0.1)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(0.1), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.1))  +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )

ggsave(filename = "Rxn with CNO.jpg",
       plot = last_plot(),
       width = 7.1, height = 5.75, units = "in", dpi = 300)  



# Frequencies ---------------------------------------------------

Rxn_Tone_On_Tone %>%
  filter(rat_ID %in% Tone_on_Tone_rats) %>%
  filter(task %in% c("Base case")) %>%
  mutate(no_go = as.factor(no_go)) %>%
  ggplot(aes(x = Position, y = Rxn_diff, shape = detail, color = no_go)) +
  # geom_point(aes(group = ID, color = Genotype), alpha = 0.3)+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(0.1), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.1))  +
  scale_color_manual(values = c("BBN" = "black", "6kHz" = "red", "12kHz" = "blue", "24kHz" = "green")) +
  facet_wrap(~ go) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )

## Individual graphs -----
individual_graphs = Rxn %>%
  filter(rat_ID %in% Tone_on_Tone_rats) %>%
  filter(task %in% c("Base case")) %>% 
  group_by(date, rat_ID, rat_name, Sex, Genotype, task, detail, phase, go, no_go) %>%
  do(mutate(., Rxn_norm = Rxn/filter(., Position == min(Position))$Rxn,
            Rxn_diff = Rxn - filter(., Position == min(Position))$Rxn)) %>%
  ungroup %>%
  mutate(no_go = as.factor(no_go)) %>%
  group_by(rat_ID) %>%
  do(
    plot = ggplot(., aes(x = Position, y = Rxn_diff, shape = detail, color = no_go)) +
      # geom_point(aes(group = ID, color = Genotype), alpha = 0.3)+
      stat_summary(fun = mean,
                   fun.min = function(x) mean(x) - se(x),
                   fun.max = function(x) mean(x) + se(x),
                   geom = "errorbar", width = 0, position = position_dodge(0.1)) +
      stat_summary(fun = mean,
                   geom = "point", position = position_dodge(0.1), size = 3) +
      stat_summary(fun = mean, geom = "line", position = position_dodge(0.1))  +
      scale_color_manual(values = c("BBN" = "black", "6kHz" = "red", "12kHz" = "blue", "24kHz" = "green")) +
      facet_wrap(~ go) +
      labs(title = glue("{unique(.$rat_name)} Tone on Tone")) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(color = "white")
      )
  )

individual_graphs$plot

# Base case by genotype ---------------------------------------------------

Rxn_table %>%
  filter(str_detect(Genotype, pattern = "Fmr1")) %>%
  filter(task %in% c("Base case", "Rotating") & detail %in% c("4-6", "Round 1")) %>%
  # filter(rat_ID < 300) %>%
  ggplot(aes(x = Position, y = Rxn_diff, color = Genotype)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 0, position = position_dodge(0.03)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(0.03), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.03), linewidth = 1)  +
  scale_x_continuous(breaks = seq(1:10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )


# Tsc2 by sex -------------------------------------------------------------

Rxn_table %>%
  filter(str_detect(Genotype, pattern = "Tsc2")) %>%
  filter(task %in% c("Base case", "Rotating") & detail %in% c("4-6", "Round 1")) %>%
  # filter(rat_ID < 300) %>%
  ggplot(aes(x = Position, y = Rxn_diff, color = Genotype, linetype = Sex)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - se(x),
               fun.max = function(x) mean(x) + se(x),
               geom = "errorbar", width = 0, position = position_dodge(0.03)) +
  stat_summary(fun = mean,
               geom = "point", position = position_dodge(0.03), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.03), linewidth = 1)  +
  scale_x_continuous(breaks = seq(1:10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )
