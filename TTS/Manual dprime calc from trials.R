
# Functions ---------------------------------------------------------------
Format_for_Psycho <- function(df) {
  check = df %>% filter(Trial_type == 0) %>% count() %>% as.numeric()
  CRnum = (if (check == 1) filter(df, Trial_type == 0) %>% .$CR %>% as.numeric() else check)
  suppressWarnings({
    FAnum = (if (check == 1) filter(df, Trial_type == 0) %>% .$FA %>% as.numeric() else check)
  })
  if(!("hit" %in% colnames(df))) df = df %>% add_column(hit = NA)
  if(!("miss" %in% colnames(df))) df = df %>% add_column(miss = NA)
  if(!("FA" %in% colnames(df))) {
    df = df %>% add_column(FA = NA)
    FAnum = 0
  }
  if(!("CR" %in% colnames(df))) {
    df = df %>% add_column(CR = NA)
    CRnum = 0
  }
  
  new_df = df %>% filter(Trial_type == 1) %>%
    mutate(CR = CRnum,
           FA = FAnum,
           hit = as.numeric(hit),
           miss = as.numeric(miss)) %>% replace(is.na(.), 0)
  return(new_df)
}

# Signal detection index calculation
Calculate_dprime <- function(df) {
  r = dprime(n_hit = df$hit,
             n_fa = df$FA,
             n_miss = df$miss,
             n_cr = df$CR,
             adjusted = TRUE)
  
   r[["bppd"]] = NULL # drop this column always (it can be wrongly dimensioned and break things)
  
  r = r %>% as_tibble() %>%
    mutate(`Inten (dB)` = df$`Inten (dB)`,
           `Freq (kHz)` = df$`Freq (kHz)`,
           `Dur (ms)` = df$`Dur (ms)`)
  
  return(r)
}



# Data prep ---------------------------------------------------------------

Daily_trial_stats = 
  core_trials %>%
  # TODO: Add filter for omit_list
  dplyr::filter(Block_number != 1) %>%
  summarise(hit = sum(Response == "Hit"),
            miss = sum(Response == "Miss"),
            CR = sum(Response == "CR"),
            FA = sum(Response == "FA"),
            .by = c(UUID, 
                    `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`, Trial_type))

Daily_trial_caluclations =
  Daily_trial_stats %>%
  nest(.by = UUID) %>%
  # .[c(1,2),] %>%
  rowwise() %>%
  mutate(formatted_data = Format_for_Psycho(data) %>% nest %>% .$data,
         dprime_table = Calculate_dprime(formatted_data) %>% nest %>% .$data)

daily_hits = select(Daily_trial_caluclations, UUID, formatted_data) %>% 
  unnest(formatted_data) 

daily_psycho = select(Daily_trial_caluclations, UUID, dprime_table) %>% 
  unnest(dprime_table)

Daily_acustic_table = 
  left_join(daily_hits,
            daily_psycho,
            by = join_by(UUID, 
                         `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`)) %>%
  left_join(core_data %>% select(date, UUID, rat_name, rat_ID, Sex, HL_state, 
                                 BG_type, BG_Intensity, file_name, phase, task, detail),
            by = join_by(UUID))


# Graph -------------------------------------------------------------------

Daily_acustic_table %>%
  filter(`Freq (kHz)` != 0) %>%
  filter(`Dur (ms)` == 50) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  # mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
  #          factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(BG_type != "White") %>%
  filter(! (BG_type == "Pink" & BG_Intensity == "30")) %>%
  mutate(Background = if_else(BG_type == "None", "None",
                              paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "White noise at 50dB", "Pink noise at 50dB"))) %>%
  mutate(`Freq (kHz)` = as.factor(`Freq (kHz)`)) %>%
  filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  ggplot(aes(x = `Inten (dB)`, y = c, linetype = HL_state, color = Background)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
  # geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
  #             # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
  # ) +
  scale_color_manual(values = c("black", "salmon")) +
  facet_wrap( ~ `Freq (kHz)`, ncol = 2) +
  labs(x = "Intensity (dB)",
       y = "beta (mean +/- SE)",
       title = "Response Bias (beta)",
       color = "Background",
       linetype = "Hearing Loss",
       shape = "Hearing Loss",
       # caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)})'"))
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  scale_y_continuous(breaks = seq(-1, 5, by = 0.5)) +
  # ylim(-1, 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )



# Functions ---------------------------------------------------------------
Format_for_Psycho <- function(df) {
  check = df %>% filter(Trial_type == 0) %>% count() %>% as.numeric()
  CRnum = (if (check == 1) filter(df, Trial_type == 0) %>% .$CR %>% as.numeric() else check)
  suppressWarnings({
    FAnum = (if (check == 1) filter(df, Trial_type == 0) %>% .$FA %>% as.numeric() else check)
  })
  if(!("hit" %in% colnames(df))) df = df %>% add_column(hit = NA)
  if(!("miss" %in% colnames(df))) df = df %>% add_column(miss = NA)
  if(!("FA" %in% colnames(df))) {
    df = df %>% add_column(FA = NA)
    FAnum = 0
  }
  if(!("CR" %in% colnames(df))) {
    df = df %>% add_column(CR = NA)
    CRnum = 0
  }
  
  new_df = df %>% filter(Trial_type == 1) %>%
    mutate(CR = CRnum,
           FA = FAnum,
           hit = as.numeric(hit),
           miss = as.numeric(miss)) %>% replace(is.na(.), 0)
  return(new_df)
}

# Signal detection index calculation
Calculate_dprime <- function(df) {
  r = dprime(n_hit = df$hit,
             n_fa = df$FA,
             n_miss = df$miss,
             n_cr = df$CR,
             adjusted = TRUE)
  
  r[["bppd"]] = NULL # drop this column always (it can be wrongly dimensioned and break things)
  
  r = r %>% as_tibble() %>%
    mutate(`Inten (dB)` = df$`Inten (dB)`,
           `Freq (kHz)` = df$`Freq (kHz)`,
           `Dur (ms)` = df$`Dur (ms)`)
  
  return(r)
}



# Data prep ---------------------------------------------------------------

Daily_trial_stats = 
  core_trials %>%
  # TODO: Add filter for omit_list
  dplyr::filter(Block_number != 1) %>%
  summarise(hit = sum(Response == "Hit"),
            miss = sum(Response == "Miss"),
            CR = sum(Response == "CR"),
            FA = sum(Response == "FA"),
            .by = c(UUID, 
                    `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`, Trial_type))

Daily_trial_caluclations =
  Daily_trial_stats %>%
  nest(.by = UUID) %>%
  # .[c(1,2),] %>%
  rowwise() %>%
  mutate(formatted_data = Format_for_Psycho(data) %>% nest %>% .$data,
         dprime_table = Calculate_dprime(formatted_data) %>% nest %>% .$data)

daily_hits = select(Daily_trial_caluclations, UUID, formatted_data) %>% 
  unnest(formatted_data) 

daily_psycho = select(Daily_trial_caluclations, UUID, dprime_table) %>% 
  unnest(dprime_table)

Daily_acustic_table = 
  left_join(daily_hits,
            daily_psycho,
            by = join_by(UUID, 
                         `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`)) %>%
  left_join(core_data %>% select(date, UUID, rat_name, rat_ID, Sex, HL_state, 
                                 BG_type, BG_Intensity, file_name, phase, task, detail),
            by = join_by(UUID))


# Graph -------------------------------------------------------------------

Daily_acustic_table %>%
  filter(`Freq (kHz)` != 0) %>%
  filter(`Dur (ms)` == 50) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  # mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
  #          factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(BG_type != "White") %>%
  filter(! (BG_type == "Pink" & BG_Intensity == "30")) %>%
  mutate(Background = if_else(BG_type == "None", "None",
                              paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "White noise at 50dB", "Pink noise at 50dB"))) %>%
  mutate(`Freq (kHz)` = as.factor(`Freq (kHz)`)) %>%
  filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  ggplot(aes(x = `Inten (dB)`, y = c, linetype = HL_state, color = Background)) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
  # geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
  #             # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
  # ) +
  scale_color_manual(values = c("black", "salmon")) +
  facet_wrap( ~ `Freq (kHz)`, ncol = 2) +
  labs(x = "Intensity (dB)",
       y = "beta (mean +/- SE)",
       title = "Response Bias (beta)",
       color = "Background",
       linetype = "Hearing Loss",
       shape = "Hearing Loss",
       # caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)})'"))
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  scale_y_continuous(breaks = seq(-1, 5, by = 0.5)) +
  # ylim(-1, 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


## dprime

# Functions ---------------------------------------------------------------
Format_for_Psycho <- function(df) {
  check = df %>% filter(Trial_type == 0) %>% count() %>% as.numeric()
  CRnum = (if (check == 1) filter(df, Trial_type == 0) %>% .$CR %>% as.numeric() else check)
  suppressWarnings({
    FAnum = (if (check == 1) filter(df, Trial_type == 0) %>% .$FA %>% as.numeric() else check)
  })
  if(!("hit" %in% colnames(df))) df = df %>% add_column(hit = NA)
  if(!("miss" %in% colnames(df))) df = df %>% add_column(miss = NA)
  if(!("FA" %in% colnames(df))) {
    df = df %>% add_column(FA = NA)
    FAnum = 0
  }
  if(!("CR" %in% colnames(df))) {
    df = df %>% add_column(CR = NA)
    CRnum = 0
  }
  
  new_df = df %>% filter(Trial_type == 1) %>%
    mutate(CR = CRnum,
           FA = FAnum,
           hit = as.numeric(hit),
           miss = as.numeric(miss)) %>% replace(is.na(.), 0)
  return(new_df)
}

# Signal detection index calculation
Calculate_dprime <- function(df) {
  r = dprime(n_hit = df$hit,
             n_fa = df$FA,
             n_miss = df$miss,
             n_cr = df$CR,
             adjusted = TRUE)
  
  r[["bppd"]] = NULL # drop this column always (it can be wrongly dimensioned and break things)
  
  r = r %>% as_tibble() %>%
    mutate(`Inten (dB)` = df$`Inten (dB)`,
           `Freq (kHz)` = df$`Freq (kHz)`,
           `Dur (ms)` = df$`Dur (ms)`)
  
  return(r)
}



# Data prep ---------------------------------------------------------------

Daily_trial_stats = 
  core_trials %>%
  # TODO: Add filter for omit_list
  dplyr::filter(Block_number != 1) %>%
  summarise(hit = sum(Response == "Hit"),
            miss = sum(Response == "Miss"),
            CR = sum(Response == "CR"),
            FA = sum(Response == "FA"),
            .by = c(UUID, 
                    `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`, Trial_type))

Daily_trial_caluclations =
  Daily_trial_stats %>%
  nest(.by = UUID) %>%
  # .[c(1,2),] %>%
  rowwise() %>%
  mutate(formatted_data = Format_for_Psycho(data) %>% nest %>% .$data,
         dprime_table = Calculate_dprime(formatted_data) %>% nest %>% .$data)

daily_hits = select(Daily_trial_caluclations, UUID, formatted_data) %>% 
  unnest(formatted_data) 

daily_psycho = select(Daily_trial_caluclations, UUID, dprime_table) %>% 
  unnest(dprime_table)

Daily_acustic_table = 
  left_join(daily_hits,
            daily_psycho,
            by = join_by(UUID, 
                         `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`)) %>%
  left_join(core_data %>% select(date, UUID, rat_name, rat_ID, Sex, HL_state, 
                                 BG_type, BG_Intensity, file_name, phase, task, detail),
            by = join_by(UUID)) 

Rat_acustic_table =
  Daily_acustic_table %>%
  reframe(hit = mean(hit, na.rm = TRUE),
          miss = mean(miss, na.rm = TRUE),
          CR = mean(CR, na.rm = TRUE),
          FA = mean(FA, na.rm = TRUE),
          dprime = mean(dprime, na.rm = TRUE),
          beta = mean(beta, na.rm = TRUE),
          aprime = mean(aprime, na.rm = TRUE),
          c = mean(c, na.rm = TRUE),
          .by = c(rat_name, rat_ID, Sex, HL_state, BG_type, BG_Intensity,
                  `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`))

Rat_acustic_table_change =
  Rat_acustic_table %>%
  filter(rat_ID %in% rats_survived_to_post_HL) %>%
    group_by(rat_name, rat_ID, Sex, BG_type, BG_Intensity,
             `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`) %>%
  do(dprime = filter(., HL_state == "baseline")$dprime - filter(., HL_state == "post-HL")$dprime)


Rat_acustic_table_change_HL =
  Rat_acustic_table %>%
  filter(rat_ID %in% rats_survived_to_post_HL) %>%
  group_by(rat_name, rat_ID, Sex, BG_type, BG_Intensity,
           `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`) %>%
  do(bind_cols(dprime = filter(., HL_state == "baseline")$dprime - filter(., HL_state == "post-HL")$dprime,
     aprime = filter(., HL_state == "baseline")$aprime - filter(., HL_state == "post-HL")$aprime,
     beta = filter(., HL_state == "baseline")$beta - filter(., HL_state == "post-HL")$beta,
     c = filter(., HL_state == "baseline")$c - filter(., HL_state == "post-HL")$c))

Rat_acustic_table_change_BG =
  Rat_acustic_table %>%
  filter(rat_ID %in% rats_survived_to_post_HL) %>%
  # get only white noise at 50dB
  filter(BG_type != "White") %>%
  filter(! (BG_type == "Pink" & BG_Intensity == "30")) %>%
  # convert BG to single variable
  mutate(Background = if_else(BG_type == "None", "None",
                              paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "Pink noise at 50dB", "White noise at 50dB"))) %>%
  # calculate difference
  group_by(rat_name, rat_ID, Sex, HL_state,
           `Freq (kHz)`, `Inten (dB)`, `Dur (ms)`) %>%
  do(bind_cols(dprime = filter(., Background == "None")$dprime - filter(., Background == "Pink noise at 50dB")$dprime,
               aprime = filter(., Background == "None")$aprime - filter(., Background == "Pink noise at 50dB")$aprime,
               beta = filter(., Background == "None")$beta - filter(., Background == "Pink noise at 50dB")$beta,
               c = filter(., Background == "None")$c - filter(., Background == "Pink noise at 50dB")$c))

# Graph -------------------------------------------------------------------

## Response bias ----
Rat_acustic_table %>%
  filter(`Freq (kHz)` != 0) %>%
  filter(`Dur (ms)` == 50) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  # mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
  #          factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(BG_type != "White") %>%
  filter(! (BG_type == "Pink" & BG_Intensity == "30")) %>%
  mutate(Background = if_else(BG_type == "None", "None",
                              paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "Pink noise at 50dB", "White noise at 50dB"))) %>%
  mutate(`Freq (kHz)` = as.factor(`Freq (kHz)`)) %>%
  filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  # beta, c or non-parametric bppd
  ggplot(aes(x = `Inten (dB)`, y = beta, linetype = HL_state, color = Background)) +
  # labels
  geom_text(label = "Convervative", color = "black", size = 4, x = 5, y = 1.8) +
  geom_text(label = "Liberal", color = "black", size = 4, x = 5, y = 0.4) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
  # geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
  #             # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
  # ) +
  scale_color_manual(values = c("black", "salmon")) +
  facet_wrap( ~ `Freq (kHz)`, ncol = 2) +
  labs(x = "Intensity (dB)",
       y = "beta (mean +/- SE)",
       title = "Response Bias (beta)",
       color = "Background",
       linetype = "Hearing Loss",
       shape = "Hearing Loss",
       # caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)})'"))
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  scale_y_continuous(breaks = seq(-1, 5, by = 0.5)) +
  # ylim(-1, 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

### response bias HL change ----
Rat_acustic_table_change_HL %>%
  filter(`Freq (kHz)` != 0) %>%
  filter(`Dur (ms)` == 50) %>%
  # filter(BG_type != "White") %>%
  # filter(! (BG_type == "Pink" & BG_Intensity == "30")) %>%
  mutate(Background = if_else(BG_type == "None", "None",
                              paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "Pink noise at 50dB", "White noise at 50dB"))) %>%
  mutate(`Freq (kHz)` = as.factor(`Freq (kHz)`)) %>%
  filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  # beta, c or non-parametric bppd
  ggplot(aes(x = `Inten (dB)`, y = beta, color = Background)) +
  # labels
  geom_text(label = "Convervative", color = "black", size = 4, x = 5, y = 0.7) +
  # geom_text(label = "Liberal", color = "black", size = 4, x = 5, y = 0.4) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 5, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
  # geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
  #             # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
  # ) +
  scale_color_manual(values = c("black", "salmon", "mediumvioletred", "bisque4")) +
  facet_wrap( ~ `Freq (kHz)`, ncol = 2) +
  labs(x = "Intensity (dB)",
       y = "beta (mean +/- SE)",
       title = "Change in Response Bias (beta) following Hearing loss",
       color = "Background",
       linetype = "Hearing Loss",
       shape = "Hearing Loss",
       # caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)})'"))
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  scale_y_continuous(breaks = seq(-1, 5, by = 0.5)) +
  # ylim(-1, 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

### Response bias BG change ----
Rat_acustic_table %>%
  filter(`Freq (kHz)` != 0) %>%
  filter(`Dur (ms)` == 50) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  mutate(HL_state = if_else(HL_state == "post-HL", "Hearing Loss", "Baseline") %>%
           factor(levels = c("Baseline", "Hearing Loss"))) %>%
  mutate(`Freq (kHz)` = as.factor(`Freq (kHz)`)) %>%
  # filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  # beta, c or non-parametric bppd
  # ggplot(aes(x = `Freq (kHz)`, y = beta, fill = HL_state)) + # By Frequency
  ggplot(aes(x = HL_state, y = beta, fill = HL_state)) + # Just HL_state
  # labels
  geom_text(label = "Convervative", color = "black", size = 4, x = 0.75, y = 2.57) +
  geom_text(label = "Liberal in background noise", color = "black", size = 4, x = 1.15, y = 0.23) +
  geom_boxplot() +
  scale_fill_manual(values = c("darkgrey", "white")) +
  labs(x = "Intensity (dB)",
       y = "beta (mean +/- SE)",
       title = "Change in Response Bias (beta) of tone-in-noise compared to quiet",
       fill = ""
       ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    # panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )
  
  
  # ggsave(filename = "beta change for HL.svg",
  #        path = save_folder,
  #        plot = last_plot(),
  #        # width = 8, height = 6, units = "in", dpi = 300)
  #        width = 5, height = 6, units = "in", dpi = 300)

## dprime ----
Rat_acustic_table %>%
  filter(`Freq (kHz)` != 0) %>%
  filter(`Dur (ms)` == 50) %>%
  filter(HL_state %in% c("baseline", "post-HL")) %>%
  # mutate(HL_state = if_else(HL_state == "post-HL", "After Hearing Loss and recovery", "Baseline") %>%
  #          factor(levels = c("Baseline", "After Hearing Loss and recovery"))) %>%
  filter(BG_type != "White") %>%
  filter(! (BG_type == "Pink" & BG_Intensity == "30")) %>%
  mutate(Background = if_else(BG_type == "None", "None",
                              paste0(BG_type, " noise at ", BG_Intensity, "dB")) %>%
           factor(levels = c("None", "Pink noise at 30dB", "Pink noise at 50dB", "White noise at 50dB"))) %>%
  mutate(`Freq (kHz)` = as.factor(`Freq (kHz)`)) %>%
  filter(! str_detect(`Inten (dB)`, pattern = "5$")) %>%
  # dprime or aprime
  ggplot(aes(x = `Inten (dB)`, y = dprime, linetype = HL_state, color = Background)) +
  geom_hline(yintercept = 1.5, linetype = "longdash", alpha = 0.5) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - FSA::se(x),
               fun.max = function(x) mean(x) + FSA::se(x),
               geom = "errorbar", width = 1, position = position_dodge(0)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0), size = 3) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0)) +
  # geom_smooth(se = FALSE, na.rm = TRUE, linewidth = 1.2,
  #             # method = "nls", formula = y ~ SSasymp(x, yf, y0, log_alpha),
  # ) +
  scale_color_manual(values = c("black", "salmon")) +
  facet_wrap( ~ `Freq (kHz)`, ncol = 2) +
  labs(x = "Intensity (dB)",
       y = "beta (mean +/- SE)",
       title = "d' (or a' is non-parametric)",
       color = "Background",
       linetype = "Hearing Loss",
       shape = "Hearing Loss",
       # caption = parse(text = glue("'Primary effect of background noise ('*chi^2*' = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$statistic, digits = 2)}, p = {round(filter(kruskal_results_dprime, str_detect(model, pattern = 'BG'))$p.adj, digits = 2)})'"))
  ) +
  scale_x_continuous(breaks = seq(-50, 90, by = 10)) +
  # ylim(-1, 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )


# Stats -------------------------------------------------------------------

## BG change ----

### data prep ----
change_in_BG.aov.data = Rat_acustic_table_change_BG %>% ungroup
change_in_BG.aov.data$beta_Gaus = LambertW::Gaussianize(change_in_BG.aov.data$beta)[, 1]
change_in_BG.aov.data$c_Gaus = LambertW::Gaussianize(change_in_BG.aov.data$beta)[, 1]


### beta aov model ----
## C is more significant than beta
  change_in_BG.bias.aov = aov(beta ~ HL_state * `Freq (kHz)`, data = change_in_BG.aov.data)
  
  Parametric_Check(change_in_BG.bias.aov)
  # Not even close to normal
  summary(change_in_BG.bias.aov)
  
  # Non-normal aov
  change_in_BG.bias.aov.Kruskal =
    tibble(
      Comparison = c(kruskal.test(beta ~ HL_state, data = change_in_BG.aov.data)$data.name,
                     kruskal.test(beta ~ `Freq (kHz)`, data = change_in_BG.aov.data)$data.name),
      P.unadj = c(kruskal.test(beta ~ HL_state, data = change_in_BG.aov.data)$p.value,
                  kruskal.test(beta ~ `Freq (kHz)`, data = change_in_BG.aov.data)$p.value)
    ) %>% 
    bind_rows(change_in_BG.aov.data %>%
                summarise(Comparison = kruskal.test(beta ~ HL_state)$data.name,
                          P.unadj = kruskal.test(beta ~ HL_state)$p.value,
                          .by = `Freq (kHz)`)) %>%
    mutate(P.adj = p.adjust(P.unadj, "BH"),
           sig = stars.pval(P.adj))
  
  change_in_BG.bias.aov.Kruskal
  
  # Post-Hoc testing
  change_in_BG.bias.aov.postHoc = 
    dunnTest(beta ~ `Freq (kHz)`, data = change_in_BG.aov.data)$res %>%
    mutate(sig = stars.pval(P.adj))
  
  filter(change_in_BG.bias.aov.postHoc, ! sig %in% c(" ", "."))


### dprime aov model ----
  change_in_BG.dprime.aov = aov(dprime ~ HL_state * `Freq (kHz)`, data = change_in_BG.aov.data)
  
  Parametric_Check(change_in_BG.dprime.aov)
  # Not even close to normal
  summary(change_in_BG.dprime.aov)
  
  # Non-normal aov
  change_in_BG.dprime.aov.Kruskal =
    tibble(
      Comparison = c(kruskal.test(dprime ~ HL_state, data = change_in_BG.aov.data)$data.name,
                     kruskal.test(dprime ~ `Freq (kHz)`, data = change_in_BG.aov.data)$data.name),
      P.unadj = c(kruskal.test(dprime ~ HL_state, data = change_in_BG.aov.data)$p.value,
                  kruskal.test(dprime ~ `Freq (kHz)`, data = change_in_BG.aov.data)$p.value)
    ) %>% 
    bind_rows(change_in_BG.aov.data %>%
                summarise(Comparison = kruskal.test(dprime ~ HL_state)$data.name,
                          P.unadj = kruskal.test(dprime ~ HL_state)$p.value,
                          .by = `Freq (kHz)`)) %>%
    mutate(P.adj = p.adjust(P.unadj, "BH"),
           sig = stars.pval(P.adj))
  
  change_in_BG.dprime.aov.Kruskal
  
  # Post-Hoc testing
  change_in_BG.dprime.aov.postHoc = 
    dunnTest(dprime ~ `Freq (kHz)`, data = change_in_BG.aov.data)$res %>%
    mutate(sig = stars.pval(P.adj))
  
  filter(change_in_BG.dprime.aov.postHoc, ! sig %in% c(" ", "."))
