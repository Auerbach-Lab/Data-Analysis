
# Load Data (optional) ----------------------------------------------------

# source("Fmr1 vs Tsc2 in LE BBN tone data.R")

# Variables ---------------------------------------------------------------

# save_folder = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Tsc2 + Rapamycin Graphs"

drop_seizure_rats = FALSE
rats_with_spontanious_seizures = c("TP3", "Blue3")

# Graphing ----------------------------------------------------------------
# Calculate standard error (SE) like standard deviation (SD)
se <- function(x, ...) {sqrt(var(x, ...)/length(x))}

n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}

# Functions ---------------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}


# N -----------------------------------------------------------------------

TH_table %>%
  filter(line == "Tsc2-LE") %>%
  select(1:5) %>%
  unique() %>%
  mutate(group = case_when(rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  summarise(n = n(), .by = c(genotype, sex, group))


# TH analysis -------------------------------------------------------------
# note that for Tukey everything needs to be a factor and Duration (a number
# value) tries to make itself continuous

#BBN only with duration & detail
Tsc2.TH.aov.data = TH_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Frequency == 0) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  filter(Duration %in% c(300, 50)) %>% # Recheck was only on 50 & 300
  mutate(genotype = as.factor(genotype),
         Duration = as.factor(Duration),
         sex = as.factor(sex))

Tsc2.TH.aov.data$Gaus = LambertW::Gaussianize(Tsc2.TH.aov.data$TH)[, 1]

Tsc2.TH.aov = aov(Gaus ~ genotype * Duration * sex,
                  data = Tsc2.TH.aov.data)

Parametric_Check(Tsc2.TH.aov)

# Non-Normal
summary(Tsc2.TH.aov)

# # Kruskal Testing - only
# lapply(c("Frequency", "genotype", "Duration"),
#        function(x) kruskal.test(reformulate(x, "TH"), data = Tsc2.TH.aov.data)) %>%
#   # Convert to table
#   do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
#   # do a p adjustment and then sig label
#   mutate(adj.p.value = p.adjust(p.value, "BH"),
#          sig = gtools::stars.pval(adj.p.value))

broom::tidy(TukeyHSD(Tsc2.TH.aov)) %>% mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  filter(sig != " " & str_detect(term, pattern = ":", negate = TRUE)) 


# No sex or genotyp effect on thresholds

# TH Averages -------------------------------------------------------------

Tsc2_TH_averages = 
  TH_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Frequency == 0) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) 

# View
Tsc2_TH_averages %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  filter(Duration %in% c(300, 50)) %>%
  group_by(genotype, Frequency, group, Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop") %>%
  # Combine Frequency and Duration to create a single key column
  mutate(key = paste0(Frequency, "kHz", "_", Duration, "ms")) %>%
  # Remove redundant columns
  select(-all_of(c("Frequency", "Duration"))) %>%
  spread(key, TH)


# TH Graph ----------------------------------------------------------------

## Threshold by group -----
Tsc2_TH_gaph = 
  TH_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Frequency == 0) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
           factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH, shape = line,
             fill = genotype, color = group, group = interaction(group, line, genotype))) +
    geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
    # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
    stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
                 position = position_dodge(1), vjust = 2, size = 3) +
    scale_color_manual(values = c("Group 1" = "darkblue", "Group 2" = "goldenrod", "Group 2 Recheck" = "green")) +
    scale_fill_manual(values = c("WT" = "grey40", "Het" = "deepskyblue", "KO" = "lightcoral")) +
    labs(x = "",
         y = "Threshold (dB, mean +/- SE)",
         caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats"),
         fill = "Genotype") +
    facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    )

print(Tsc2_TH_gaph)

TH_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Frequency == 0) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
           factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH, shape = line,
             fill = genotype, color = group, group = interaction(group, line, genotype))) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  # geom_point(aes(color = genotype), alpha = 0.3, position = position_dodge(1)) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
               position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("Group 1" = "darkblue", "Group 2" = "goldenrod", "Group 2 Recheck" = "green")) +
  scale_fill_manual(values = c("WT" = "grey40", "Het" = "deepskyblue", "KO" = "lightcoral")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats"),
       fill = "Genotype") +
  facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Threshold by Frequency -----
TH_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Duration == 50) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  #filter(group == "Group 1") %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
           factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH, fill = genotype, group = genotype)) +
  geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
  scale_color_manual(values = c("Group 1" = "darkblue", "Group 2" = "goldenrod", "Group 2 Recheck" = "green")) +
  scale_fill_manual(values = c("WT" = "grey40", "Het" = "deepskyblue", "KO" = "lightcoral")) +
  labs(x = "",
       y = "Threshold (dB, mean +/- SE)",
       fill = "Genotype") +
  facet_wrap( ~ Frequency, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = c(.2,.2)
  )

# ggsave(filename = "Tsc2_TH_all_freq.svg",
#        path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Tsc2 + Rapamycin Graphs",
#        plot = last_plot(),
#        width = 8, height = 8, units = "in", dpi = 150)


## Overall Threshold ----
Tsc2_TH_overview_gaph = 
  TH_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Frequency == 0) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>% 
           factor(levels = c("BBN", "4", "8", "16", "32"))) %>%
  ggplot(aes(x = genotype, y = TH, fill = genotype, group = interaction(genotype))) +
    geom_boxplot(position = position_dodge(1), linewidth = 1, width = 0.8) +
    stat_summary(aes(color = sex, group = interaction(sex, genotype)),
                 fun = function(x) mean(x, na.rm = TRUE),
                 geom = "point", show.legend = FALSE, size = 3) +
    stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, 
                 position = position_dodge(1), vjust = 2, size = 3) +
    # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
    scale_color_manual(values = c("Male" = "blue", "Female" = "coral")) +
    scale_fill_manual(values = c("WT" = "darkgrey", "Het" = "deepskyblue", "KO" = "lightcoral")) +
    labs(x = "",
         y = "Threshold (dB, mean +/- SE)",
         caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats"),
         fill = "Tsc2\nGenotype") +
    facet_wrap( ~ Duration, ncol = 5, scales = "free_x") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    )

print(Tsc2_TH_overview_gaph)

## Threshold by group/sex ---
TH_table %>%
  filter(line == "Tsc2-LE") %>%
  filter(Frequency == 0) %>%
  filter(! detail %in% c("Mixed", "Rotating")) %>%
  filter(Duration %in% c(300)) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2 bad original",
                           .default = "Unknown"),
         Duration = factor(Duration, levels = c(50, 100, 300), ordered = TRUE)) %>%
  filter(group %in% c("Group 1", "Group 2")) %>%
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
  facet_wrap(~ sex) +
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

# Rxn analysis ------------------------------------------------------------
## Get data ----
Tsc2_Rxn = Rxn_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck"))

Tsc2_Rxn$Gaus = LambertW::Gaussianize(Tsc2_Rxn$Rxn)[, 1]


## BBN Model ----
  Tsc2.Rxn.BBN.aov.data = Tsc2_Rxn %>%
    filter(Frequency == 0) %>%
    # filter(sex == "Male") %>%
    filter(Duration %in% c(50, 100, 300)) %>%
    filter(! str_detect(Intensity, pattern = "5$")) %>% # group 1 didn't have 5 steps at 10dB diff
    filter(Intensity < 90 & Intensity > 10) # corrects for measuring differences between groups
    
  Tsc2.Rxn.BBN.aov = aov(Gaus ~ sex * Duration * genotype,
                         data = Tsc2.Rxn.BBN.aov.data)
  
  # Normality testing
  Parametric_Check(Tsc2.Rxn.BBN.aov)
  
  # Not even close to normal
  # summary(Tsc2.Rxn.BBN.aov)
  
  # Kruskal Testing - Main effects only 
  lapply(c("genotype", "sex", "Duration", # Main effects
           "interaction(genotype, sex)", 
           "interaction(genotype, Duration)", 
           "interaction(sex, Duration)",
           "interaction(genotype, sex, Duration)"
         ), 
         function(x) kruskal.test(reformulate(x, "Rxn"),
                                  data = filter(Tsc2.Rxn.BBN.aov.data))) %>% 
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "bonf"),
           sig = gtools::stars.pval(adj.p.value)) %>%
    select(method, parameter, statistic, data.name, adj.p.value, sig)

  # Post-Hoc Dunn's Test
  FSA::dunnTest(Rxn ~ interaction(Duration, genotype, sex),
                  data = Tsc2.Rxn.BBN.aov.data,
                  method = "bonf") %>%
    .$res %>%
    as_tibble() %>%
    select(-P.unadj) %>%
    mutate(Sig = gtools::stars.pval(P.adj),
           Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
           Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
           dur1 = str_split_fixed(Comp1, '\\.', 3)[,1] %>% as.numeric(),
           geno1 = str_split_fixed(Comp1, '\\.', 3)[,2],
           sex1 = str_split_fixed(Comp1, '\\.', 3)[,3],
           dur2 = str_split_fixed(Comp2, '\\.', 3)[,1] %>% as.numeric(),
           geno2 = str_split_fixed(Comp2, '\\.', 3)[,2],
           sex2 = str_split_fixed(Comp2, '\\.', 3)[,3]) %>%
    # only compare within a sex (sib-sib direct comparison)
    # filter(sex1 == sex2) %>%
    filter(! Sig %in% c(" ", ".")) %>%
    arrange(dur1) %>%
    select(-Comparison, -Comp1, -Comp2)
  
## 50ms only -----
  Tsc2.Rxn.BBN50.aov.data = Tsc2_Rxn %>%
    filter(Frequency == 0) %>%
    filter(Duration %in% c(50)) %>%
    filter(! str_detect(Intensity, pattern = "5$")) %>% # group 1 didn't have 5 steps at 10dB diff
    filter(Intensity < 90 & Intensity > 10) # corrects for measuring differences between groups
  
  Tsc2.Rxn.BBN50.aov = aov(Gaus ~ sex * genotype,
                            data = Tsc2.Rxn.BBN300.aov.data)
  
  # Normality testing
  Parametric_Check(Tsc2.Rxn.BBN50.aov)
  
  # Not even close to normal
  # summary(Tsc2.Rxn.BBN.aov)
  
  # Kruskal Testing - Main effects only 
  lapply(c("genotype", "sex", # Main effects
           "interaction(genotype, sex)"
  ), 
  function(x) kruskal.test(reformulate(x, "Rxn"),
                           data = filter(Tsc2.Rxn.BBN50.aov.data))) %>% 
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "bonf"),
           sig = gtools::stars.pval(adj.p.value)) %>%
    select(method, parameter, statistic, data.name, adj.p.value, sig)
  
  # Post-Hoc Dunn's Test
  FSA::dunnTest(Rxn ~ interaction(genotype, sex),
                data = Tsc2.Rxn.BBN50.aov.data,
                method = "bonf") %>%
    .$res %>%
    as_tibble() %>%
    select(-P.unadj) %>%
    mutate(Sig = gtools::stars.pval(P.adj),
           Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
           Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
           geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
           sex1 = str_split_fixed(Comp1, '\\.', 2)[,2],
           geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
           sex2 = str_split_fixed(Comp2, '\\.', 2)[,2]) %>%
    # only compare within a sex (sib-sib direct comparison)
    # filter(sex1 == sex2) %>%
    filter(! Sig %in% c(" ", ".")) %>%
    select(-Comparison, -Comp1, -Comp2)

## 300ms only -----
  Tsc2.Rxn.BBN300.aov.data = Tsc2_Rxn %>%
    filter(Frequency == 0) %>%
    filter(Duration %in% c(300)) %>%
    filter(! str_detect(Intensity, pattern = "5$")) %>% # group 1 didn't have 5 steps at 10dB diff
    filter(Intensity < 90 & Intensity > 10) # corrects for measuring differences between groups
  
  Tsc2.Rxn.300BBN.aov = aov(Gaus ~ sex * genotype,
                         data = Tsc2.Rxn.BBN300.aov.data)
  
  # Normality testing
  Parametric_Check(Tsc2.Rxn.300BBN.aov)
  
  # Not even close to normal
  # summary(Tsc2.Rxn.BBN.aov)
  
  # Kruskal Testing - Main effects only 
  lapply(c("genotype", "sex", # Main effects
           "interaction(genotype, sex)"
  ), 
  function(x) kruskal.test(reformulate(x, "Rxn"),
                           data = filter(Tsc2.Rxn.BBN300.aov.data))) %>% 
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "bonf"),
           sig = gtools::stars.pval(adj.p.value)) %>%
    select(method, parameter, statistic, data.name, adj.p.value, sig)
  
  # Post-Hoc Dunn's Test
  FSA::dunnTest(Rxn ~ interaction(genotype, sex),
                data = Tsc2.Rxn.BBN300.aov.data,
                method = "bonf") %>%
    .$res %>%
    as_tibble() %>%
    select(-P.unadj) %>%
    mutate(Sig = gtools::stars.pval(P.adj),
           Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
           Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
           geno1 = str_split_fixed(Comp1, '\\.', 2)[,1],
           sex1 = str_split_fixed(Comp1, '\\.', 2)[,2],
           geno2 = str_split_fixed(Comp2, '\\.', 2)[,1],
           sex2 = str_split_fixed(Comp2, '\\.', 2)[,2]) %>%
    # only compare within a sex (sib-sib direct comparison)
    # filter(sex1 == sex2) %>%
    filter(! Sig %in% c(" ", ".")) %>%
    select(-Comparison, -Comp1, -Comp2)

# Rxn Graphs --------------------------------------------------------------

## Durations on BBN ----
Tsc_single_frequency_rxn_graph =   
  Rxn_table %>%
    filter(line == "Tsc2-LE") %>%
    {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
    filter(detail %in% c("Alone", "Recheck")) %>%
    mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                             rat_ID < 300 ~ "Group 1",
                             rat_ID >= 300 ~ "Group 2",
                             .default = "Unknown")) %>%
    filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
    filter(Duration %in% c(300, 50)) %>%
    filter(Frequency == 0) %>%
    mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
    filter(! str_detect(Intensity, pattern = "5$")) %>%
    filter(Intensity < 90 & Intensity > 10) %>%
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
         color = "Genotype", linetype = "",
         caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats")) +
    scale_linetype_manual(values = c("Group 1" = "solid", "Group 2 Recheck" = "longdash")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
    scale_x_continuous(breaks = seq(0, 90, by = 10)) +
    facet_wrap(~ Duration) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    ) 

print(Tsc_single_frequency_rxn_graph)

# ggsave(filename = "Tsc2_Rxn_all_freq.jpg",
# path = save_folder,
# plot = last_plot(),
# width = 5, height = 6, units = "in", dpi = 300)

# TH_annotation = 

## Tones ----
TH_annotation = 
  TH_table %>%
    filter(line == "Tsc2-LE") %>%
    {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
    filter(Frequency != 0) %>%
    filter(detail %in% c("Alone", "Recheck")) %>%
    mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                             rat_ID < 300 ~ "Group 1",
                             rat_ID >= 300 ~ "Group 2",
                             .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  group_by(genotype, Frequency) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")


Rxn_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  filter(sex == "Male") %>%
  filter(Frequency != 0) %>%
  filter(Intensity <= 90 & Intensity >= 10) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = genotype, group = genotype)) +
  ## Lines for each group
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 1.5, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line", position = position_dodge(1)) +
  # threshold lines
  geom_vline(data = TH_annotation, 
             aes(xintercept = TH, color = genotype, group = genotype), 
             linetype = "longdash", show.legend = FALSE) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "",
       caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats")) +
  scale_linetype_manual(values = c("Group 1" = "solid", "Group 2 Recheck" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Frequency) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

ggsave(filename = "Tsc2_Rxn_tones_MALES.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Tsc2 + Rapamycin Graphs",
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 150)


## Sex differences graph ---------------------------------------------------

Tsc_sex_differences_rxn_graph =   
  Rxn_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  filter(Duration %in% c(300, 50)) %>%
  filter(Frequency == 0) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == "BBN") %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(sex),
             color = genotype, group = interaction(sex, genotype))) +
  ## Overall average lines
  # stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
  #              fun = function(x) mean(x, na.rm = TRUE),
  #              fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
  #              fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
  #              geom = "errorbar", width = 1.5, linetype = "solid", linewidth = 1.5, alpha = 0.5,
  #              position = position_dodge(4)) +
  # stat_summary(aes(x = Intensity, y = Rxn,color = genotype,group = genotype),
  #              fun = function(x) mean(x, na.rm = TRUE),
  #              geom = "line", linetype = "solid", linewidth = 1.5, alpha = 0.5,
  #              position = position_dodge(4)) +
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
       color = "Genotype", linetype = "",
       caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats")) +
  # scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

print(Tsc_sex_differences_rxn_graph)

### Sex diff: individual females on graph ----

Rxn_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  filter(Duration %in% c(300, 50)) %>%
  filter(sex == "Female") %>%
  filter(Frequency == 0) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == "BBN") %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(sex),
             color = genotype, group = interaction(sex, genotype))) +
  ## Individuals
  # geom_line(aes(group = rat_ID), alpha = 0.5) +
  # geom_dl(aes(label = rat_name), method = list(dl.combine("last.points"), cex = 0.8),
  #         hjust = 3) +
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
       color = "Genotype", linetype = "",
       caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats")) +
  # scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

### Sex diff: individual males on graph ----

Rxn_table %>%
  filter(line == "Tsc2-LE") %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  filter(Duration %in% c(300, 50)) %>%
  filter(sex == "Male") %>%
  filter(Frequency == 0) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN")) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  filter(Intensity < 90 & Intensity > 10) %>%
  filter(Frequency == "BBN") %>%
  ggplot(aes(x = Intensity, y = Rxn, linetype = as.factor(sex),
             color = genotype, group = interaction(sex, genotype))) +
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
       color = "Genotype", linetype = "",
       caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats")) +
  # scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

# Individual Graphs -------------------------------------------------------

BBN_Individual_Graphs = 
Rxn_table %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else .} %>%
  filter(line == "Tsc2-LE") %>%
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

# d' ----------------------------------------------------------------------

## Get data ----
dprime_table_Tsc2 = core_data %>%
  # Omit Training & Reset days
  filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  filter(line == "Tsc2-LE") %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                           rat_ID < 300 ~ "Group 1",
                           rat_ID >= 300 ~ "Group 2",
                           .default = "Unknown")) %>%
  filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
  # Get dprimes
  unnest(dprime) %>%
  summarise(dprime = mean(dprime, na.rm = TRUE), 
            .by = c(rat_name, rat_ID, genotype, sex, group, Freq, dB, Dur)) %>%
  rename(Frequency = Freq, Intensity = dB, Duration = Dur)

## BBN d' analysis ----
Tsc2.dprime.BBN.aov.data = dprime_table_Tsc2 %>%
  filter(Frequency == 0) %>%
  # filter(sex == "Male") %>%
  filter(Duration %in% c(50, 100, 300)) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>% # group 1 didn't have 5 steps at 10dB diff
  filter(Intensity < 90 & Intensity > 10) # corrects for measuring differences between groups


Tsc2.dprime.BBN.aov.data$Gaus = LambertW::Gaussianize(Tsc2.dprime.BBN.aov.data$dprime)[, 1]


Tsc2.dprime.BBN.aov = aov(Gaus ~ sex * Duration * genotype,
                          data = Tsc2.dprime.BBN.aov.data)

# Normality testing
Parametric_Check(Tsc2.dprime.BBN.aov)

# Not even close to normal

# Kruskal Testing - Main effects only 
lapply(c("genotype", "sex", "Duration", # Main effects
         "interaction(genotype, sex)", 
         "interaction(genotype, Duration)", 
         "interaction(sex, Duration)",
         "interaction(genotype, sex, Duration)"
), 
function(x) kruskal.test(reformulate(x, "dprime"),
                         data = filter(Tsc2.dprime.BBN.aov.data))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, adj.p.value, sig)

# Post-Hoc Dunn's Test
FSA::dunnTest(dprime ~ interaction(Duration, genotype, sex),
              data = Tsc2.dprime.BBN.aov.data,
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         dur1 = str_split_fixed(Comp1, '\\.', 3)[,1] %>% as.numeric(),
         geno1 = str_split_fixed(Comp1, '\\.', 3)[,2],
         sex1 = str_split_fixed(Comp1, '\\.', 3)[,3],
         dur2 = str_split_fixed(Comp2, '\\.', 3)[,1] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, '\\.', 3)[,2],
         sex2 = str_split_fixed(Comp2, '\\.', 3)[,3]) %>%
  # only compare within a sex (sib-sib direct comparison)
  # filter(sex1 == sex2) %>%
  filter(! Sig %in% c(" ", ".")) %>%
  arrange(dur1) %>%
  select(-Comparison, -Comp1, -Comp2)

## Graph ----
dprime_table_Tsc2 %>%
  {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
  filter(Frequency == 0) %>%
  filter(Duration != 100) %>%
  filter(Intensity <= 50 & Intensity >= 10) %>%
  filter(! str_detect(Intensity, pattern = "5$")) %>%
  ggplot(aes(x = Intensity, y = dprime, color = genotype, group = genotype)) +
  ## Lines for each group
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - se(x),
               fun.max = function(x) mean(x, na.rm = TRUE) + se(x),
               geom = "errorbar", width = 1.5, position = position_dodge(1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(1), size = 3) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "line", position = position_dodge(1)) +
  # threshold lines
  geom_vline(data = TH_table %>%
               filter(line == "Tsc2-LE") %>%
               {if (drop_seizure_rats) filter(., !(rat_name %in% rats_with_spontanious_seizures)) else (.)} %>%
               filter(Frequency == 0) %>%
               filter(Duration != 100) %>%
               filter(detail %in% c("Alone", "Recheck")) %>%
               mutate(group = case_when(detail == "Recheck" ~ "Group 2 Recheck",
                                        rat_ID < 300 ~ "Group 1",
                                        rat_ID >= 300 ~ "Group 2",
                                        .default = "Unknown")) %>%
               filter(group %in% c("Group 1", "Group 2 Recheck")) %>%
               group_by(genotype, Duration, sex) %>%
               summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop"), 
             aes(xintercept = TH, color = genotype, group = genotype), 
             linetype = "longdash", show.legend = FALSE) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "Genotype", linetype = "",
       caption = if_else(drop_seizure_rats, "Without Female rats with spontanious seizures", "All rats")) +
  scale_linetype_manual(values = c("Group 1" = "solid", "Group 2 Recheck" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  facet_wrap(Duration ~ sex) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 


ggsave(filename = "Tsc2_dprime_BBN.svg",
       path = save_folder,
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 150)
