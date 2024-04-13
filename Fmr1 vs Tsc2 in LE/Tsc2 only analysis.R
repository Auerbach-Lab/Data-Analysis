
# Load Data (optional) ----------------------------------------------------

# source("Fmr1 vs Tsc2 in LE BBN tone data.R")

# Variables ---------------------------------------------------------------

drop_TP3 = TRUE

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


# TH Averages -------------------------------------------------------------
Tsc2_TH_averages = TH_table %>%
  filter(line == "Tsc2-LE") %>%
  group_by(genotype, detail, Frequency , Duration) %>%
  summarise(TH = mean(TH, na.rm = TRUE), .groups = "drop")


# # TH analysis -------------------------------------------------------------
# # note that for Tukey everything needs to be a factor and Duration (a number
# # value) tries to make itself continuous
# 
# Tsc2.TH.aov.data = TH_table %>%
#   filter(line == "Tsc2-LE") %>%
#   filter(Duration == "50" & detail == "Alone") %>%
#   # Getting an NA from TP5 so drop
#   filter(! (rat_name == "TP5" & Frequency == "4"))
# 
# Tsc2.TH.aov.data$Gaus = LambertW::Gaussianize(Tsc2.TH.aov.data$TH)[, 1]
# 
# Tsc2.TH.aov = aov(Gaus ~ genotype * as.factor(Frequency), 
#                   data = Tsc2.TH.aov.data)
# 
# Parametric_Check(Tsc2.TH.aov)
# 
# # Non-Normal
# # summary(Tsc2.TH.aov)
# 
# # Kruskal Testing - only 
# lapply(c("Frequency", "genotype"), 
#        function(x) kruskal.test(reformulate(x, "TH"), data = Tsc2.TH.aov.data)) %>% 
#   # Convert to table
#   do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
#   # do a p adjustment and then sig label
#   mutate(adj.p.value = p.adjust(p.value, "BH"),
#          sig = gtools::stars.pval(adj.p.value))
# 
# # No effect on thresholds
# 
# # TODO: BBN only with duration & detail
# 
# Tsc2_TH_averages %>%
#   filter(Frequency == 0) %>%
#   # Combine Frequency and Duration to create a single key column
#   mutate(key = paste0(Frequency, "kHz", "_", Duration, "ms")) %>%
#   # Remove redundant columns
#   select(-all_of(c("Frequency", "Duration"))) %>%
#   spread(key, TH)
# 
# # TH Duration --------------------------------------------------------------
# 
# Tsc2.TH.BBN.aov.data = TH_table %>%
#   filter(line == "Tsc2-LE") %>%
#   filter(detail != "Rotating" & Frequency == 0) 
# 
# Tsc2.TH.BBN.aov.data$Gaus = LambertW::Gaussianize(Tsc2.TH.BBN.aov.data$TH)[, 1]
# 
# Tsc2.TH.BBN.aov = aov(TH ~ detail * Duration * genotype,
#                       data = Tsc2.TH.BBN.aov.data)
# 
# Parametric_Check(Tsc2.TH.BBN.aov)
# 
# # Normal
# summary(Tsc2.TH.BBN.aov)
# 
# # Only primary effects so no post-Hoc needed


# TH Graph ----------------------------------------------------------------

Tsc2_TH_gaph = 
  TH_table %>%
    filter(line == "Tsc2-LE") %>%
    filter(Frequency == 0) %>%
    filter(detail %in% c("Alone", "Recheck")) %>%
    {if (drop_TP3) filter(., rat_name != "TP3")} %>%
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
    # scale_color_manual(values = c("Tsc2-LE" = "darkblue", "Fmr1-LE" = "red")) +
    scale_color_manual(values = c("Group 1" = "darkblue", "Group 2" = "goldenrod", "Group 2 Recheck" = "green")) +
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

print(Tsc2_TH_gaph)

# Rxn analysis ------------------------------------------------------------
# Get data
Tsc2_Rxn_over_TH = filter(Rxn_table_over_TH, line == "Tsc2-LE")

Tsc2_Rxn_over_TH$Gaus = LambertW::Gaussianize(Tsc2_Rxn_over_TH$Rxn)[, 1]


## Overall Model
  Tsc2_Rxn_overall_model = aov(Gaus ~ Frequency * Intensity * genotype * sex,
                               data = filter(Tsc2_Rxn_over_TH, detail == "Alone" & Duration == "50"))
  # Normality testing
  Parametric_Check(Tsc2_Rxn_overall_model)
  
  # Non-normal even with transformation
  kruskal.test(Rxn ~ genotype,
               data = filter(Tsc2_Rxn_over_TH, detail == "Alone" & Duration == "50"))
  
  kruskal.test(Rxn ~ Frequency,
               data = filter(Tsc2_Rxn_over_TH, detail == "Alone" & Duration == "50"))


## BBN Model
  Tsc2.Rxn.BBN.aov.data = Tsc2_Rxn_over_TH %>%
    filter(Frequency == 0)
  
  Tsc2.Rxn.BBN.aov = aov(Gaus ~ sex * Duration * genotype,
                         data = Tsc2.Rxn.BBN.aov.data)
  
  # Normality testing
  Parametric_Check(Tsc2.Rxn.BBN.aov)
  
  # Not even close to normal
  # summary(Tsc2.Rxn.BBN.aov)
  
  # Kruskal Testing - only 
  lapply(c("genotype", "sex", "interaction(genotype, sex)"), 
         function(x) kruskal.test(reformulate(x, "Rxn"),
                                  data = filter(Tsc2.Rxn.BBN.aov.data, detail == "Alone"))) %>% 
    # Convert to table
    do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
    # do a p adjustment and then sig label
    mutate(adj.p.value = p.adjust(p.value, "bonf"),
           sig = gtools::stars.pval(adj.p.value))

  # Post-Hoc Dunn's Test
  FSA::dunnTest(Rxn ~ interaction(Duration, genotype, sex), 
                  data = filter(Tsc2.Rxn.BBN.aov.data, detail == "Alone"),
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
    filter(dur1 == dur2) %>%
    filter(! Sig %in% c(" ", ".")) %>%
    arrange(dur1) %>%
    select(-Comparison, -Comp1, -Comp2, -dur2) 
    
  



# Rxn Graphs --------------------------------------------------------------

single_Frequency = "BBN"

Tsc_single_frequency_rxn_graph =   
  Rxn_table %>%
    {if (drop_TP3) filter(., rat_name != "TP3")} %>%
    filter(line == "Tsc2-LE") %>%
    filter(Duration %in% c(50, 100, 300)) %>%
    # rename(Intensity = `Inten (dB)`) %>%
    filter(detail %in% c("Alone", "Recheck")) %>%
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
         color = "Genotype", linetype = "",
         caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
    scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
    scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
    scale_x_continuous(breaks = seq(0, 90, by = 10)) +
    facet_wrap(~ Duration) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
    ) 

print(Tsc_single_frequency_rxn_graph)

# ggsave(filename = paste0("Tsc2_Rxn_all_freq.jpg"),
# plot = last_plot(),
# width = 5, height = 6, units = "in", dpi = 300)


# Sex differences graph ---------------------------------------------------

Tsc_sex_differences_rxn_graph =   
  Rxn_table %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
  filter(line == "Tsc2-LE") %>%
  filter(Duration %in% c(50, 100, 300)) %>%
  # rename(Intensity = `Inten (dB)`) %>%
  filter(detail %in% c("Alone", "Recheck")) %>%
  mutate(group = if_else(rat_ID < 300, "Group 1", "Group 2")) %>%
  # filter(rat_ID < 314) %>%
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
       caption = if_else(drop_TP3, "Without Het F TP3", "With TP3")) +
  # scale_linetype_manual(values = c("Group 1" = "solid", "Group 2" = "longdash")) +
  scale_color_manual(values = c("WT" = "black", "Het" = "deepskyblue", "KO" = "red")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  # facet_wrap(~ Duration) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

print(Tsc_sex_differences_rxn_graph)


# Individual Graphs -------------------------------------------------------

BBN_Individual_Graphs = 
Rxn_table %>%
  {if (drop_TP3) filter(., rat_name != "TP3")} %>%
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
