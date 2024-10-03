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
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
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
    plot = ggplot(., aes(x = Position, y = Rxn_diff*1000, shape = detail, color = no_go)) +
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
      ylim(-100, 50) +
      scale_x_continuous(breaks = seq(1, 10, by = 1)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(color = "white")
      )
  )

individual_graphs$plot

# Rxn for Base case by genotype -------------------------------------------

Rxn_table %>%
  filter(str_detect(Genotype, pattern = "Fmr1")) %>%
  filter(task %in% c("Base case", "Rotating") & detail %in% c("4-6", "Round 1")) %>%
  mutate(group = case_when(rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
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


# Probe histogram ---------------------------------------------------------

Probe_table %>%
  filter(rat_ID > 300) %>%
  filter(Response == "FA") %>%
  filter(str_detect(Genotype, pattern ="Fmr1-LE")) %>%
  mutate(Genotype = str_remove(Genotype, pattern = "Fmr1-LE_")) %>%
  ggplot(aes(x = `Reaction_(s)`, y = Genotype, color = Genotype, group = rat_ID)) +
  geom_rect(aes(xmin = 2.25, xmax = 3, ymin = -Inf, ymax = Inf), 
            fill = "lightgrey", color = "lightgrey", show.legend = FALSE) + # position 6
  geom_rect(aes(xmin = 3, xmax = 3.75, ymin = -Inf, ymax = Inf), 
            fill = "darkgrey", color = "darkgrey", show.legend = FALSE) + # position 5
  geom_rect(aes(xmin = 3.75, xmax = 4.5, ymin = -Inf, ymax = Inf), 
            fill = "lightgrey", color = "lightgrey", show.legend = FALSE) + # position 6
  geom_vline(xintercept = c(0.05, 0.8, 1.55, 2.3, 3.05, 3.8, 4.55, 5.3), color = "grey50") + # tone ends
  #geom_vline(xintercept = c(0, 0.75, 1.5, 2.25, 3, 3.75, 4.5, 5.25), color = "grey50") + # response window starts
  geom_beeswarm(cex = 2, method = "compactswarm") +
  scale_x_continuous(breaks = seq(from = 0, to = 6, by = 1), expand = c(0.01, 0)) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  labs(x = "Reaction time (s)",
       y = "Genotype (Fmr1-LE rats)",
       color = "Genotype") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )

## FAs on non-probe trials ----
Probe_FA_count %>%
  filter(rat_ID > 300) %>%
  ggplot(aes(x = position, y = FA_count, color = Genotype)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 0, position = position_dodge(0.03)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(0.03), size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line", position = position_dodge(0.03), linewidth = 1)  +
  scale_x_continuous(breaks = seq(1:10)) +
  theme_classic() +
  theme(
    # legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "white")
  )


# FA by genotype ----------------------------------------------------------


