# Functions ---------------------------------------------------------------
n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}


# Get Data ----------------------------------------------------------------

# Get pilot BBN threshold data
BBN_alone_TH_pilot = 
  TH_table_detail %>%
  # limit to individuals with baseline & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # remove group 3 (Blues) that have no temporal integration & in-progress group
  filter(rat_ID < 200) %>%
  # limit to BBN with no background noise
  filter(detail == "Alone" & Frequency == 0 & BG_Intensity  == "None")

# Get pilot BBN reaction time data
BBN_alone_Rxn_pilot =
  Rxn_table_detail %>%
  # limit to individuals with baseine & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # remove group 3 (Blues) that have no temporal integration & in-progress group
  filter(rat_ID < 200) %>%
  # limit to BBN with no background noise
  filter(detail == "Alone" & Frequency == 0 & BG_Intensity  == "None") 

# All data
BBN_TempInt_TH = 
  TH_table_detail %>%
  # limit to individuals with baseline & post-HL
  filter(! rat_name %in% c("Green2", "Orange1", "Blue2")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # # Remove struggling rat
  # filter(rat_name != "Green24") %>%
  # limit to BBN with no background noise
  filter(detail %in% c("Alone", "Recheck", "Recheck 2") & Frequency == 0 & BG_Intensity  == "None") %>%
  # Add group numbers to look for effects between groups
  mutate(group = case_when(rat_ID %in% Group_TTS_pilot ~ "Pilot",
                           rat_ID %in% Group_1 ~ "Group 1",
                           rat_ID %in% Group_2 ~ "Group 2",
                           rat_ID %in% Group_3 & detail == "Recheck 2" ~ "Group 3 Recheck 2",
                           rat_ID %in% Group_3 & detail == "Recheck" ~ "Group 3 Redo",
                           rat_ID %in% Group_3 ~ "Group 3",
                           .default = "Unknown") %>% ordered(levels = c("Pilot", "Group 1", "Group 2", "Group 3", "Group 3 Redo", "Group 3 Recheck 2")))

BBN_TempInt_Rxn =
  Rxn_table_detail %>%
  # limit to individuals with baseine & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # limit to BBN with no background noise
  filter(detail %in% c("Alone", "Recheck", "Recheck 2") & Frequency == 0 & BG_Intensity  == "None") %>%
  # Add group numbers to look for effects between groups
  mutate(group = case_when(rat_ID %in% Group_TTS_pilot ~ "Pilot",
                           rat_ID %in% Group_1 ~ "Group 1",
                           rat_ID %in% Group_2 ~ "Group 2",
                           rat_ID %in% Group_3 & detail == "Recheck 2" ~ "Group 3 Recheck 2",
                           rat_ID %in% Group_3 & detail == "Recheck" ~ "Group 3 Redo",
                           rat_ID %in% Group_3 ~ "Group 3",
                           .default = "Unknown") %>% ordered(levels = c("Pilot", "Group 1", "Group 2", "Group 3", "Group 3 Redo", "Group 3 Recheck 2")))

BBN_TempInt_Rxn_Mixed = 
  Rxn_table_detail %>%
  # limit to individuals with baseine & post-HL
  filter(! rat_name %in% c("Green2", "Orange1")) %>%
  # Remove rat with permanent threshold shifts
  filter(rat_ID != 191 | rat_name != "Green12") %>%
  # limit to BBN with no background noise
  filter(detail %in% c("Mixed") & Frequency == 0 & BG_Intensity  == "None") %>%
  # Add group numbers to look for effects between groups
  mutate(group = case_when(rat_ID %in% Group_TTS_pilot ~ "Pilot",
                           rat_ID %in% Group_1 ~ "Group 1",
                           rat_ID %in% Group_2 ~ "Group 2",
                           rat_ID %in% Group_3 & detail == "Recheck 2" ~ "Group 3 Recheck 2",
                           rat_ID %in% Group_3 & detail == "Recheck" ~ "Group 3 Redo",
                           rat_ID %in% Group_3 ~ "Group 3",
                           .default = "Unknown") %>% ordered(levels = c("Pilot", "Group 1", "Group 2", "Group 3", "Group 3 Redo", "Group 3 Recheck 2")))


# Calculate Change in Rxn -------------------------------------------------

BBN_TempInt_Rxn_change =
  BBN_TempInt_Rxn %>%
  group_by(rat_ID, rat_name, group, Frequency, Duration, Intensity) %>%
  do(Rxn = filter(., HL_state == "post-HL")$Rxn - filter(., HL_state == "baseline")$Rxn %>%
       ifelse(identical(., numeric(0)), NA_integer_, .)) %>% 
  mutate(Rxn = ifelse(identical(Rxn, numeric(0)), NA_integer_, Rxn))


# Export for Ben ----------------------------------------------------------
# 
# BBN_alone_TH_pilot %>% 
#     fwrite(paste0(save_folder, "BBN_alone_TH_pilot_", Sys.Date(),".csv"), row.names = FALSE)
# 
# BBN_alone_Rxn_pilot %>% 
#   fwrite(paste0(save_folder, "BBN_alone_Rxn_pilot_", Sys.Date(),".csv"), row.names = FALSE)


# Recheck Thresholds ------------------------------------------------------

Threshold_Recheck = 
  BBN_TempInt_TH %>% 
    filter(rat_ID %in% Group_3) %>%
    group_by(rat_ID, rat_name, Duration) %>%
    do(tibble(TH_1 = filter(., ! detail %in% c("Recheck", "Recheck 2"))$TH, 
              TH_2 = if(! is_na(filter(., detail == "Recheck")$TH)) filter(., detail == "Recheck")$TH
                        else NA_integer_,
              TH_3 = if(! is_na(filter(., detail == "Recheck 2")$TH)) filter(., detail == "Recheck 2")$TH
                        else NA_integer_,
              Improvement = round(TH_1 - TH_3, digits = 1))
       ) %>%
  arrange(Duration, rat_ID)

print(Threshold_Recheck)

# Threshold Graph ---------------------------------------------------------

# # Boxplot - Pilot group pre and post noise exposure
# BBN_alone_TH_pilot %>%
# # BBN_TempInt_TH %>%
#   mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
#            factor(levels = c("4", "8", "16", "32", "BBN")),
#          BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
#          HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
#   filter(HL_state %in% c("baseline", "post-HL")) %>%
#   mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
#            factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
#   ggplot(aes(x = as.factor(Duration), y = TH, fill = HL_state, group = interaction(Duration, HL_state))) +
#     geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
#     stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
#     scale_fill_manual(values = c("#F8766D", "#00BFC4"), guide = "legend") +
#     labs(title = "BBN thresholds by durations",
#          x = "Stimulus Duration (ms)",
#          y = "Threshold (dB)",
#          fill = "Condition",
#          caption = "Dropped individual with permanent threshold shifts. We only have data for the 50ms alone for all rats."
#          ) +
#     # facet_wrap( ~ HL_state, ncol = 5, scales = "free_x") +
#     theme_classic() +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       legend.position = c(0.85, 0.93),
#       legend.background=element_blank(),
#       panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#     )

# ggsave(filename = "TH_BBN_pilot.jpg",
#        path = save_folder,
#        plot = last_plot(),
#        width = 8, height = 6, units = "in", dpi = 300)

# BBN_TempInt_TH %>%
#   mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
#            factor(levels = c("4", "8", "16", "32", "BBN")),
#          BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
#          HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
#   filter(HL_state %in% c("baseline", "post-HL")) %>%
#   mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
#            factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
#   ggplot(aes(x = as.factor(Duration), y = TH,
#              color = HL_state, fill = group,
#              group = interaction(Duration, group, HL_state))) +
#   geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
#   stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
#   scale_color_manual(values = c("coral", "grey30"), guide = "legend") +
#   labs(title = "BBN thresholds by durations",
#        x = "Stimulus Duration (ms)",
#        y = "Threshold (dB)"
#   ) +
#   # facet_wrap( ~ HL_state, ncol = 5, scales = "free_x") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     # legend.position = c(0.85, 0.93),
#     legend.background=element_blank(),
#     panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#   )

TempInt_TH_graph = 
BBN_TempInt_TH %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
           factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(group %in% c("Group 1", "Group 3", "Group 3 Redo", "Group 3 Recheck 2")) %>%
  filter(HL_state == "Baseline") %>%
  ggplot(aes(x = as.factor(Duration), y = TH,
             color = HL_state, fill = group,
             group = interaction(Duration, group, HL_state))) +
  geom_boxplot(na.rm = TRUE, position = position_dodge(1), linewidth = 1, width = 0.8) +
  stat_summary(fun.data = n_fun, geom = "text", show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 2) +
  scale_color_manual(values = c("coral", "grey30"), guide = "legend") +
  labs(title = "BBN thresholds by durations",
       x = "Stimulus Duration (ms)",
       y = "Threshold (dB)",
       color = "Hearing loss?"
  ) +
  # facet_wrap( ~ HL_state, ncol = 5, scales = "free_x") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.position = c(0.85, 0.93),
    legend.background=element_blank(),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  )

print(TempInt_TH_graph)

# Reaction Graph ----------------------------------------------------------

BBN_alone_Rxn_pilot %>%
# BBN_TempInt_Rxn %>%
  #limit to individuals with all 3 durations
  # filter(rat_name %in% c("Orange11", "Orange12", "Green11")) %>%
  filter(Intensity < 90) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  ggplot(aes(x = Intensity, y = Rxn, color = HL_state, group = interaction(HL_state))) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
  geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
              # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha)
  ) +
  # geom_text(data = Rxn_annotations, aes(label = sig), size = 10, show.legend = FALSE) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       color = "",
       title = "Reaction curves for BBN before and after Hearing Loss",
       caption = "Following hearing loss, there is a loss of temproal integration. We only have data for the 50ms alone for all rats."
       ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  ) +
  facet_wrap( ~ Duration, ncol = 3, scales = "fixed") +
  theme(legend.position = c(0.9, 0.8),
        legend.background=element_blank())

# ggsave(filename = glue("BBN_alone_Rxn_pilot.jpg"),
#        path = save_folder,
#        plot = last_plot(),
#        width = 8, height = 6, units = "in", dpi = 300)


TempInt_Rxn_Graph =
BBN_TempInt_Rxn %>%
  filter(Intensity < 95 & Intensity > 19) %>%
  # filter(rat_name != "Green21") %>%
  filter(HL_state == "baseline") %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN")),
         BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
         HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  # filter(group %in% c("Group 1", "Group 3", "Group 2", "Group 3 Redo")) %>%
  filter(group %in% c("Group 3", "Group 3 Redo", "Group 3 Recheck 2")) %>%
  ggplot(aes(x = Intensity, y = Rxn, 
             color = as.factor(Duration), shape = as.factor(Duration),
             group = interaction(HL_state, as.factor(Duration)))) +
    # Individual lines
    # geom_line(aes(linetype = as.factor(rat_name),
    #               group = interaction(HL_state, as.factor(Duration), rat_ID)),
    #           alpha = 0.5) +
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - FSA::se(x),
                 fun.max = function(x) mean(x) + FSA::se(x),
                 geom = "errorbar", width = 1, position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
    geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
    labs(x = "Intensity (dB)",
         y = "Reaction time (ms, mean +/- SE)",
         color = "Duration", shape = "Duration",
         title = "Reaction curves for BBN showing temporal integration at baseline",
         linetype = "Rat"
         ) +
    scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
      panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    ) +
    facet_wrap( ~ group, ncol = 2, scales = "fixed") +
    theme(#legend.position = c(0.9, 0.8),
          legend.background=element_blank())

print(TempInt_Rxn_Graph)

# BBN_TempInt_Rxn %>%
#   filter(Intensity < 90 & Intensity > 19) %>%
#   filter(HL_state == "post-HL") %>%
#   mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
#            factor(levels = c("4", "8", "16", "32", "BBN")),
#          BG = if_else(BG_type == "None", "None", paste0(BG_type, "\n", BG_Intensity, "dB")),
#          HL_state = factor(HL_state, levels = c("baseline", "HL", "recovery", "post-HL"))) %>%
#   filter(HL_state %in% c("baseline", "post-HL")) %>%
#   ggplot(aes(x = Intensity, y = Rxn,
#              color = as.factor(Duration), shape = as.factor(Duration),
#              group = interaction(HL_state, as.factor(Duration)))) +
#   # # Individual lines
#   # geom_line(aes(group = interaction(HL_state, as.factor(Duration), rat_ID))) +
#   stat_summary(fun = mean,
#                fun.min = function(x) mean(x) - FSA::se(x),
#                fun.max = function(x) mean(x) + FSA::se(x),
#                geom = "errorbar", width = 1, position = position_dodge(1)) +
#   stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
#   geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
#   labs(x = "Intensity (dB)",
#        y = "Reaction time (ms, mean +/- SE)",
#        color = "Duration", shape = "Duration",
#        title = "Reaction curves for BBN for temporal integration after Hearing Loss"
#   ) +
#   scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#     panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#   ) +
#   facet_wrap( ~ group, ncol = 2, scales = "fixed") +
#   theme(legend.position = c(0.9, 0.8),
#         legend.background=element_blank())

# BBN_TempInt_Rxn_change %>%
#   filter(Intensity < 90 & Intensity > 29) %>%
#   mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
#            factor(levels = c("4", "8", "16", "32", "BBN"))) %>%
#   ggplot(aes(x = Intensity, y = Rxn,
#              color = as.factor(Duration), shape = as.factor(Duration),
#              group = interaction(as.factor(Duration)))) +
#   stat_summary(fun = mean,
#                fun.min = function(x) mean(x) - FSA::se(x),
#                fun.max = function(x) mean(x) + FSA::se(x),
#                geom = "errorbar", width = 1, position = position_dodge(1)) +
#   stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
#   geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
#   labs(x = "Intensity (dB)",
#        y = "Change in Reaction time (ms, mean +/- SE)",
#        color = "Duration", shape = "Duration",
#        title = "Increase in reaction time after Hearing Loss"
#   ) +
#   scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#     panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
#   ) +
#   facet_wrap( ~ group, ncol = 2, scales = "fixed") +
#   theme(#legend.position = c(0.9, 0.8),
#         legend.background=element_blank())

# Individual Graphs -------------------------------------------------------

BBN_TempInt_Rxn_individual_graphs =
  BBN_TempInt_Rxn %>%
  filter(Intensity < 85 & Intensity > 19) %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN"))) %>%
  filter(group %in% c("Group 2", "Group 3", "Group 3 Redo", "Group 3 Recheck 2")) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  group_by(rat_ID, rat_name) %>%
  do(
    single_rat_graph =
    ggplot(data = .,
           aes(x = Intensity, y = Rxn, 
               color = as.factor(Duration), 
               linetype = group, shape = group,
               group = interaction(group, as.factor(Duration)))) +
    # Individual lines
    stat_summary(fun = mean,
                 fun.min = function(x) mean(x) - FSA::se(x),
                 fun.max = function(x) mean(x) + FSA::se(x),
                 geom = "errorbar", width = 1, position = position_dodge(1)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
    geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1) +
    labs(x = "Intensity (dB)",
         y = "Reaction time (ms, mean +/- SE)",
         color = "Duration", 
         shape = "Group", linetype = "Group",
         title = glue("{unique(.$rat_name)} - Reaction curves for BBN at baseline"),
    ) +
    scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
    facet_wrap( ~ HL_state, ncol = 2, scales = "fixed") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
      panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    ) +
    theme(#legend.position = c(0.9, 0.8),
      legend.background=element_blank())
  )

# Show all the individual graphs
BBN_TempInt_Rxn_individual_graphs$single_rat_graph

## Mixed Duration ------------------------------------------------------

BBN_TempInt_Rxn_Mixed_individual_graphs =
  BBN_TempInt_Rxn_Mixed %>%
  mutate(Frequency = str_replace_all(Frequency, pattern = "0", replacement = "BBN") %>%
           factor(levels = c("4", "8", "16", "32", "BBN"))) %>%
  filter(group %in% c("Group 3", "Group 3 Redo")) %>%
  filter(HL_state %in% c("baseline")) %>%
  group_by(rat_ID, rat_name) %>%
  do(
    single_rat_graph =
      ggplot(data = .,
             aes(x = Intensity, y = Rxn, 
                 color = as.factor(Duration), 
                 linetype = group, shape = group,
                 group = interaction(group, as.factor(Duration)))) +
      # Individual lines
      stat_summary(fun = mean,
                   fun.min = function(x) mean(x) - FSA::se(x),
                   fun.max = function(x) mean(x) + FSA::se(x),
                   geom = "errorbar", width = 1, position = position_dodge(1)) +
      stat_summary(fun = mean, geom = "point", position = position_dodge(1), size = 2) +
      stat_summary(fun = mean, geom = "line", position = position_dodge(1)) +
      labs(x = "Intensity (dB)",
           y = "Reaction time (ms, mean +/- SE)",
           color = "Duration", 
           shape = "Group", linetype = "Group",
           title = glue("{unique(.$rat_name)} - Reaction curves for BBN at baseline"),
      ) +
      scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
      facet_wrap( ~ HL_state, ncol = 2, scales = "fixed") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
        panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
      ) +
      theme(#legend.position = c(0.9, 0.8),
        legend.background=element_blank())
  )

# Show all the individual graphs
BBN_TempInt_Rxn_Mixed_individual_graphs$single_rat_graph


# Temp. Integration Testing -----------------------------------------------

## data import -----
load(glue("{projects_folder}/run_archive.Rdata"), .GlobalEnv)

TI_testing_core_data = 
  run_archive %>%
    # Omit Invalid runs
    filter(invalid != "TRUE") %>%
    #Omit runs with wrong delay window, the negate means it returns non-matches
    #ISSUE: gives warning because it expects a vector not a tibble
    suppressWarnings(
      filter(str_detect(as.vector(warnings_list), pattern = "wrong delay window", negate = TRUE)) 
    ) %>%
    unnest_wider(assignment) %>%
    # Only keep relevant Experiments
    filter(experiment %in% c("TTS")) %>%
    filter(rat_ID > 300) %>%
    filter(task %in% c("Duration Testing", "Rxn", "TH") & detail == "Mixed") %>%
    mutate(data_set = case_when(task == "Duration Testing" ~ "Single Intensity Testing",
                                task %in% c("Rxn", "TH") ~ "Mixed Intensity BBN")) %>%
    # record date of hearing loss, Sex, & Genotype
    left_join(rat_archive, by = c("rat_ID" = "Rat_ID")) %>%
    unnest_wider(stats) %>%
    # Omit days with > 45% FA, i.e. guessing
    filter(FA_percent < FA_cutoff) 

TI_testing_data = TI_testing_core_data %>%
  select(-reaction) %>%
  unnest(dprime) %>%
  left_join(TI_testing_core_data %>% select(-dprime) %>% unnest(reaction),
            # For Reaction & dprime
            by = join_by(!!!intersect(TI_testing_core_data %>% select(-reaction) %>% unnest(dprime) %>% names,
                                      TI_testing_core_data %>% select(-dprime) %>% unnest(reaction) %>% names),
                         Freq == `Freq (kHz)`, dB == `Inten (dB)`, Dur == `Dur (ms)`)) %>%
  reframe(Reaction = mean(Rxn * 1000, na.rm = TRUE),
          dprime = mean(dprime, na.rm = TRUE),
          .by = c(rat_name, rat_ID, data_set, Freq, dB, Dur)) %>%
  group_by(rat_name, rat_ID, Freq, dB) %>%
  do(mutate(., Rxn_norm = Reaction/filter(., Dur == 50)$Reaction,
            Rxn_diff = Reaction - filter(., Dur == 50)$Reaction))

## Single intensity Testing Graph -------------------------------------------

TI_testing_data %>%
  filter(data_set == "Single Intensity Testing") %>%
  mutate(dB = as.factor(dB)) %>%
  ggplot(aes(x = Dur, y = Rxn_diff, fill = dB, group = interaction(dB, Dur))) +
    geom_boxplot() +
    # geom_line(aes(linetype = rat_name, color = dB, group = interaction(rat_name, dB)),
    #           linewidth = 0.8) +
    theme_ipsum_es() +
    theme(legend.position = "bottom")

## Individual By duration Graph --------------------------------------------

TI_testing_data %>%
  filter(data_set == "Single Intensity Testing") %>%
  mutate(dB = as.factor(dB)) %>%
  ggplot(aes(x = Dur, y = Rxn_diff, fill = dB, group = interaction(dB, Dur))) +
    geom_line(aes(linetype = rat_name, color = dB, group = interaction(rat_name, dB)),
              linewidth = 0.8) +
    geom_dl(aes(label = rat_name), method = list("last.points","bumpup", cex = 0.8)) +
    xlim(NA, 320) +
    facet_wrap(~ dB) +
    theme_ipsum_es() +
    theme(legend.position = "bottom")

## Mixed Graphs ----------------------------------------------------------

### Rxn Graph ------
TI_testing_core_data %>%
  # filter(data_set == "Mixed Intensity BBN") %>%
  unnest(reaction) %>%
  ggplot(aes(x = `Inten (dB)`, y = Rxn, 
             color = as.factor(`Dur (ms)`), shape = as.factor(`Dur (ms)`),
             group = as.factor(`Dur (ms)`), rat_ID)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 2, position = position_dodge(0.1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.1), size = 2) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.1)) +
  labs(x = "Intensity (dB)",
       y = "Reaction time (ms, mean +/- SE)",
       title= "Reaction Time",
       color = "Duration", shape = "Duration",
       linetype = "Rat"
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10), limits = c(20, NA)) +
  facet_wrap(~ rat_name, scales = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  ) +
  theme(#legend.position = c(0.9, 0.8),
    legend.background=element_blank())

### d' Graph ------
TI_testing_core_data %>%
  # filter(data_set == "Mixed Intensity BBN") %>%
  unnest(dprime) %>%
  ggplot(aes(x = dB, y = dprime, 
             color = as.factor(Dur), shape = as.factor(Dur),
             group = as.factor(Dur), rat_ID)) +
  geom_hline(yintercept = 1.5, color = "blue", linetype = "longdash") + 
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 2, position = position_dodge(0.1)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.1), size = 2) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.1)) +
  labs(x = "Intensity (dB)",
       y = "d' (mean +/- SE)",
       title= "d'",
       color = "Duration", shape = "Duration",
       linetype = "Rat"
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  xlim(NA, 60) + 
  facet_wrap(~ rat_name, scales = "free_y") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    panel.grid.minor.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
  ) +
  theme(#legend.position = c(0.9, 0.8),
    legend.background=element_blank())

