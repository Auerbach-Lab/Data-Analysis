
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
  ggplot(aes(x = `Inten (dB)`, y = beta, linetype = HL_state, color = Background)) +
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

