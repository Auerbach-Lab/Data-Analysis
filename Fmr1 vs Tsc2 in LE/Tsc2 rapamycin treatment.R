# Load Data (optional) ----------------------------------------------------

source("Tsc2 rapamycin data.R")

# Vehicle data checks ----------------------------------------------------

## TH table ----
TH_table %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(Duration == 50) %>%
  filter(detail %in% c("Recheck", "Vehicle (Tween 80)", "None", "Post Vehicle")) %>%
  mutate(detail = str_replace(detail, ("Recheck|None"), "Baseline")) %>%
  reframe(TH = mean(TH, na.rm = TRUE),
          n = length(unique(rat_ID)),
          IDs = str_flatten_comma(unique(rat_name)),
          .by = c(genotype, detail)) %>%
  arrange(detail, genotype)


## Rxn Graph ====
Rapa_rxn_data_limited %>%
  filter(! rat_name %in% c("Lime1", "Lime2")) %>%
  filter(detail %in% c("Baseline", "Vehicle")) %>%
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
  scale_linetype_manual(values = c("Baseline" = "solid", "Vehicle" = "dotted")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    legend.position = c(.8,.85),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

### stats ----
Vehicle_rapa_aov = 
  aov(Rxn ~ genotype * detail * sex,
      data = Rapa_rxn_data_limited %>%
        filter(! rat_name %in% c("Lime1", "Lime2")) %>%
        filter(detail %in% c("Baseline", "Vehicle")))

shapiro.test(Vehicle_rapa_aov$residuals)$p.value

summary(Vehicle_rapa_aov)


### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("`Inten (dB)`", "genotype", "detail", "sex" # Main effects
), 
function(x) kruskal.test(reformulate(x, "Rxn"),
                         data = Rapa_rxn_data_limited %>%
                           filter(! rat_name %in% c("Lime1", "Lime2")) %>%
                           filter(detail %in% c("Baseline", "Vehicle")))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)


# Rapamycin treatment 1 -----------------------------------------------------
## Threshold -----

### data ----
TH_rapa_data = 
  TH_table %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
  filter(Duration == 50) %>%
  filter(rat_name != "Lime3") %>%
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)",
                       "Recheck", "Post Treatment", "None", "3+w Post Treatment")) %>%
  mutate(detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
         detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)",
                                            "Rapamycin (6mg/kg)", 
                                            "Post Treatment", "3+w Post Treatment"),
                         labels = c("Baseline", "Vehicle", "Rapamycin",
                                    "Recovery", "Permanent")))

TH_rapa_data %>%
  reframe(TH = mean(TH, na.rm = TRUE),
          .by = c(genotype, detail)) %>%
  arrange(detail, genotype)

### stats ----
TH_rapa_aov = 
  aov(TH ~ genotype * detail * sex,
      data = TH_rapa_data)

shapiro.test(TH_rapa_aov$residuals)$p.value

summary(TH_rapa_aov)

### TH Genotype * Sex ----
aov(TH ~ genotype * sex,
       data = TH_table %>%
      filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
      filter(Duration == 50) %>%
      filter(detail %in% c("Vehicle (Tween 80)",
                           "Recheck", "None")) %>%
      mutate(detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
             detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)",
                                                "Rapamycin (6mg/kg)", 
                                                "Post Treatment", "3+w Post Treatment"),
                             labels = c("Baseline", "Vehicle", "Rapamycin",
                                        "Recovery", "Permanent")))) %>%
  TukeyHSD() %>%
  tidy %>%
  mutate(sig = gtools::stars.pval(adj.p.value))

### TH Genotype * Sex Graph ----
ggplot(data = TH_table %>%
         filter(rat_ID %in% Tsc2_rapamycin_treated_rats) %>%
         filter(Duration == 50) %>%
         filter(detail %in% c("Vehicle (Tween 80)",
                              "Recheck", "None")) %>%
         mutate(detail = str_replace(detail, pattern = "None", replacement = "Recheck"),
                detail = factor(detail, levels = c("Recheck", "Vehicle (Tween 80)",
                                                   "Rapamycin (6mg/kg)", 
                                                   "Post Treatment", "3+w Post Treatment"),
                                labels = c("Baseline", "Vehicle", "Rapamycin",
                                           "Recovery", "Permanent"))),
       aes(x = genotype, y = TH, fill = genotype, color = sex)) +
  geom_boxplot(linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE,
               position = position_dodge(1), vjust = 2, size = 3) +
  labs(x = "",
       y = "Threshold (dB)",
       fill = "Genotype",
       color = "Sex") +
  scale_y_continuous(limits = c(19, 31), breaks = c(seq(18, 32, 2))) +
  scale_fill_manual(values = c("WT" = "darkgrey", "Het" = "deepskyblue", "KO" = "red")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
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


### Threshold Graph -----
TH_graph = 
  ggplot(data = TH_rapa_data,
         aes(x = genotype, y = TH, fill = genotype)) +
  geom_boxplot(linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE,
               position = position_dodge(1), vjust = 2, size = 3) +
  labs(x = "",
       y = "Threshold (dB)",
       fill = "Genotype",
       color = "Sex") +
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

print(TH_graph)

ggsave(filename = "TH_rapa.jpg",
       plot = TH_graph,
       width = 900, height = 900, units = "px", dpi = 100)

## Rxn ----
### stats ----

Rxn_rapa_aov = 
  aov(Rxn ~ `Inten (dB)` * genotype * detail * sex,
      data = Rapa_rxn_data_limited %>%
        filter(rat_name != "Lime3") %>%
        filter(detail %in% c( "Baseline", "Vehicle", "Rapamycin", "Recovery", "Permanent")))

shapiro.test(Rxn_rapa_aov$residuals)$p.value

summary(Rxn_rapa_aov)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("`Inten (dB)`", "genotype", "detail", "sex" # Main effects
), 
function(x) kruskal.test(reformulate(x, "Rxn"),
                         data = Rapa_rxn_data_limited %>%
                           filter(rat_name != "Lime3") %>%
                           filter(detail %in% c( "Baseline", "Vehicle", 
                                                 "Rapamycin", "Recovery", "Permanent")))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)


### Graph (Rapa) -----
Rapa_rxn_data_limited %>%
  filter(rat_name != "Lime3") %>%
  filter(detail %in% c("Rapamycin", "Baseline")) %>%
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
  scale_linetype_manual(values = c("Baseline" = "solid", "Vehicle " = "dotted", 
                                   "Rapamycin" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 

### Graph - Recovery -----
Rapa_rxn_data_limited %>%
  filter(rat_name != "Lime3") %>%
  filter(detail %in% c("Recovery", "Baseline")) %>%
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
  scale_linetype_manual(values = c("Baseline" = "solid", "Recovery" = "dotdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 

### Graph - Permanent -----
Rapa_rxn_data_limited %>%
  filter(rat_name != "Lime3") %>%
  filter(detail %in% c("Permanent", "Baseline")) %>%
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
  scale_linetype_manual(values = c("Baseline" = "solid", "Permanent" = "twodash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 

### Graph - Rapa Change -----
Rapa_rxn_data_limited %>%
  filter(! rat_name %in% c("Lime3")) %>%
  filter(`Inten (dB)` > 15) %>%
  filter(detail %in% c("Baseline", "Rapamycin")) %>%
  group_by(rat_ID, rat_name, sex, genotype, line, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  do(
    mutate(., Rxn_length = length(Rxn),
           Baseline = ifelse(length(Rxn) == 2, 
                             filter(., detail == "Baseline")$Rxn, NA_integer_)) %>% #print %>%
      mutate(Rxn_diff = ifelse(Rxn != Baseline,
                               (Rxn - Baseline),
                               NA_integer_)) #%>% print
  ) %>%
  ungroup %>%
  # filter(! is.na(Rxn_diff)) %>%
  filter(detail == "Rapamycin") %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn_diff, linetype = as.factor(detail),
             color = genotype, group = interaction(detail, genotype))) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 1, position = position_dodge(0.5)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(0.5)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", position = position_dodge(0.5)) +
  geom_hline(aes(yintercept = 0)) + 
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "Treatment") +
  scale_linetype_manual(values = c("Rapamycin" = "longdash", "Recovery" = "dotdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 

### Graph - Recovery Change -----
Rapa_rxn_data_limited %>%
  filter(! rat_name %in% c("Lime3")) %>%
  filter(`Inten (dB)` > 15) %>%
  filter(detail %in% c("Baseline", "Recovery")) %>%
  group_by(rat_ID, rat_name, sex, genotype, line, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  do(
    mutate(., Rxn_length = length(Rxn),
           Baseline = ifelse(length(Rxn) == 2, 
                             filter(., detail == "Baseline")$Rxn, NA_integer_)) %>% #print %>%
      mutate(Rxn_diff = ifelse(Rxn != Baseline,
                               (Rxn - Baseline),
                               NA_integer_)) #%>% print
  ) %>%
  ungroup %>%
  # filter(! is.na(Rxn_diff)) %>%
  filter(detail == "Recovery") %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn_diff, linetype = as.factor(detail),
             color = genotype, group = interaction(detail, genotype))) +
  geom_hline(aes(yintercept = 0)) + 
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
  scale_linetype_manual(values = c("Rapamycin" = "longdash", "Recovery" = "dotdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  )


### Graph - Permanent Change -----
Rapa_rxn_data_limited %>%
  filter(! rat_name %in% c("Lime3")) %>%
  filter(`Inten (dB)` > 15) %>%
  filter(detail %in% c("Baseline", "Permanent")) %>%
  group_by(rat_ID, rat_name, sex, genotype, line, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  do(
    mutate(., Rxn_length = length(Rxn),
           Baseline = ifelse(length(Rxn) == 2, 
                             filter(., detail == "Baseline")$Rxn, NA_integer_)) %>% #print %>%
      mutate(Rxn_diff = ifelse(Rxn != Baseline,
                               (Rxn - Baseline),
                               NA_integer_)) #%>% print
  ) %>%
  ungroup %>%
  # filter(! is.na(Rxn_diff)) %>%
  filter(detail == "Permanent") %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn_diff, linetype = as.factor(detail),
             color = genotype, group = interaction(detail, genotype))) +
  geom_hline(aes(yintercept = 0)) + 
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
  scale_linetype_manual(values = c("Rapamycin" = "longdash", "Permanent" = "twodash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 


# WT graph ----------------------------------------------------------------
Rapa_rxn_data_limited %>%
  filter(genotype == "WT") %>%
  filter(detail %in% c( "Baseline", "Vehicle", "Rapamycin", "Recovery", "Permanent")) %>%
  filter(`Inten (dB)` > 15) %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn, linetype = as.factor(detail),
             color = as.factor(detail), group = interaction(detail, genotype))) +
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
       color = "Treatment", linetype = "Treatment") +
  scale_color_manual(values = c("Baseline" = "black", "Vehicle" = "blue",
                                "Rapamycin" = "red", "Recovery" = "darkviolet",
                                "Permanent" = "darkgreen")) + 
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 


# Het graph ---------------------------------------------------------------
Rapa_rxn_data_limited%>%
  filter(genotype == "Het") %>%
  filter(detail %in% c( "Baseline", "Vehicle", "Rapamycin", "Recovery", "Permanent")) %>%
  filter(`Inten (dB)` > 15) %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn, linetype = as.factor(detail),
             color = as.factor(detail), group = interaction(detail, genotype))) +
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
       color = "Treatment", linetype = "Treatment") +
  scale_color_manual(values = c("Baseline" = "black", "Vehicle" = "blue",
                                "Rapamycin" = "red", "Recovery" = "darkviolet",
                                "Permanent" = "darkgreen")) + 
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.9,.85)
  ) 

# Individual Graphs -------------------------------------------------------
## Group 1----
Individual_Graphs_Gp1 = 
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

Individual_Graphs_Gp1[c(1:8),]$single_rat_graph


## Group 2 ----
Individual_Graphs_Gp2 = 
  core_data %>%
  filter(rat_ID %in% Tsc2_rapamycin_treated_rats_group2) %>%
  filter(! task %in% c("Training", "Reset")) %>%    # Omit Training & Reset days
  filter(detail %in% c("Rapamycin (6mg/kg)", "Vehicle (Tween 80)", "Post Vehicle",
                       "Post Treatment", "3+w Post Treatment", "None")) %>%
  # filter(FA_percent < FA_cutoff) %>%    # Omit days with > 45% FA, i.e. guessing
  unnest(reaction) %>% 
  filter(`Inten (dB)` >= 15) %>%
  mutate(name = rat_name,
         detail = factor(detail, ordered = TRUE,
                         levels = c("None", "Vehicle (Tween 80)", "Post Vehicle", 
                                    "Rapamycin (6mg/kg)", 
                                    "Post Treatment", "3+w Post Treatment"),
                         labels = c("Pre-Treatment", "Vehicle", "Post Vehicle",
                                    "Rapamycin", "Recovery", "Post Treatment"))) %>%
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
                                     "Vehicle" = "cyan", 
                                     "Post Vehicle" = "darkblue",
                                     "Rapamycin" = "red",
                                     "Recovery" = "goldenrod",
                                     "Post Treatment" = "forestgreen")) +
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

Individual_Graphs_Gp2[c(8, 1:4),]$single_rat_graph
