# Load Necessary Datasets -------------------------------------------------
load(glue("{projects_folder}/run_archive.Rdata"), .GlobalEnv)
rat_archive = fread(glue("{projects_folder}/rat_archive.csv"), 
                    select = c("Rat_ID", "DOB", "Sex", "Genotype", "HL_date"))

# Individual Trial Data
trials_data = fread(paste0(projects_folder, "Fmr1 SD_archive.csv.gz")) %>%
  as.tibble()



# Get core data -----------------------------------------------------------
cat("filtering...")


core_columns = c("date", "rat_name", "rat_ID",
                 "file_name", "experiment", "phase", "task", "detail",
                 "stim_type", "analysis_type", "block_size", "complete_block_count",
                 "dprime", "reaction", "FA_percent", "trial_count", "hit_percent")

#TODO: deal with multiple runs in a day, 
dataset = run_archive %>%
  # Omit Invalid runs, if any marked but they shouldn't be
  filter(invalid != "TRUE") %>%
  #Omit runs with wrong delay window, the negate means it returns non-matches
  #ISSUE: gives warning because it expects a vector not a tibble
  suppressWarnings(
    filter(str_detect(as.vector(warnings_list), pattern = "wrong delay window", negate = TRUE)) 
  ) %>%
  # Get essential columns in usable form; expands the dataframe
  unnest_wider(assignment) %>%
  # only keep the Fmr1 SD experimental data for the Octave phase
  filter(experiment == "Fmr1 SD" & phase == "Octave") %>% 
  unnest_wider(stats) %>%
  unnest(dprime) %>%
  # Decode genotype & record date of hearing loss
  left_join(rat_archive, by = c("rat_ID" = "Rat_ID")) %>%
  # Calculate time since HL
  mutate(date = ymd(date), HL_date = ymd(HL_date),
         HL = case_when(is.na(HL_date) ~ -1,
                        date <= HL_date ~ -1,
                        date > HL_date ~ as.numeric(date - HL_date),
                        TRUE ~ -100),
         HL_state = case_when(HL == -1 ~ "baseline",
                              HL < 3 ~ "HL",
                              HL >= 3 & HL < 15  ~ "recovery",
                              HL >= 15 ~ "post-HL",
                              HL == -100 ~ "issue",
                              TRUE ~ "issue") %>% factor(levels = c("baseline", "HL", "recovery", "post-HL"), ordered = TRUE)) %>%
  # Change Genotype to drop line
  mutate(Genotype = str_remove(Genotype, "Fmr1_SD_"))

Discrimination_data = dataset %>%
  # Omit Training & Reset days; this is the same as filtering to task == Discrimination
  filter(analysis_type == "Octave") %>%
  # pointless column
  select(-dprime) %>%
  # get reaction time
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  # get d'
  unnest(FA_detailed) %>%
  # rename Frequency
  rename(Frequency = `Freq (kHz)`) %>%
  # calculate octave steps - string extract w/ regex is to get the go frequency
  mutate(Go_freq = as.numeric(str_extract(file_name, pattern = "[:digit:]+?(?=-.+?kHz)")),
         NoGo_freq = str_extract(file_name, pattern = "[:digit:]+?(kHz)") %>% str_remove("kHz") %>% as.numeric(),
         octave_fraction = log(as.numeric(str_extract(file_name, pattern = "[:digit:]+?(?=-.+?kHz)"))/Frequency)/log(2),
         octave_steps = abs(round(octave_fraction * 12))) %>%
  # determine if 1/8 scale or zoomed in 1/16 scale
  mutate(Type = case_when(all((octave_steps %% 2) == 0) ~ "Broad",
                          any((octave_steps %% 2) != 0) ~ "Zoom",
                          .default = "Error"),
         Range = R.utils::seqToHumanReadable(octave_steps) %>% str_extract(pattern = "[:digit:]+-[:digit:]+"),
         .by = c(date, rat_ID)) %>%
  filter(Type != "Error")

Training_data = dataset %>%
  # Only keep Training & Holding days
  filter(analysis_type != "Octave") %>%
  filter(task == "Training") %>%
  filter(HL_state == "baseline")

#TODO: stopped here


# Write out for Walker ----------------------------------------------------

# Discrimination_data %>% select(-summary, -threshold, -comments, -warnings_list) %>%
#   fwrite(glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Walker/Fmr1_SD_Octave_", str_remove_all(Sys.Date(), "-"),".csv"), row.names = FALSE)

# Graphing ----------------------------------------------------------------
n_fun <- function(x){
  # print(x)
  return(data.frame(y = min(x), label = paste0("n = ", length(x))))
}


# Parametric test ---------------------------------------------------------
Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}

# Trial counts ------------------------------------------------------------

## Daily Trail Count by Freq --- 
# trial_daily_count = 
trials_data %>%
  filter(UUID %in% (Discrimination_data$UUID)) %>%
  left_join(Discrimination_data %>%
              select(date, rat_ID, rat_name, Genotype, HL_state,
                     task, detail, UUID, omit_list,
                     Frequency, octave_fraction, octave_steps) %>%
              unique(),
            by = join_by("UUID", `Freq (kHz)` == Frequency)) %>%
  # filter out omitted trials
  group_by(UUID) %>% 
  do(filter(., ! Trial_number %in% .$omit_list)) %>%
  ungroup %>%
  # Apply filters
  filter(detail == "Normal") %>%
  filter(HL_state == "baseline") %>%
  filter(! is.na(octave_steps)) %>%
  # get daily 
  reframe(trial_count = n(),
          Trial_type = unique(Trial_type),
          .by = c(rat_ID, rat_name, Genotype, task, detail, 
                  `Freq (kHz)`, octave_fraction, octave_steps)) %>%
  # # get rat averages 
  # reframe(avg_trial_count = mean(trial_count, na.rm = TRUE),
  #         .by = c(rat_ID, rat_name, Genotype, task, detail, 
  #                 `Freq (kHz)`, octave_fraction, octave_steps)) %>%
  # Get Averages
  summarise(avg_trial_count = mean(trial_count, na.rm = TRUE),
            SD = sd(trial_count, na.rm = TRUE),
            SE = FSA::se(trial_count, na.rm = TRUE),
            .by = c(Genotype, task, detail)) %>%
  arrange(Genotype)

### Stats -----
t.test(trial_count ~ Genotype,
       data = trials_data %>%
         filter(UUID %in% (Discrimination_data$UUID)) %>%
         left_join(Discrimination_data %>%
                     select(date, rat_ID, rat_name, Genotype, HL_state,
                            task, detail, UUID, omit_list,
                            Frequency, octave_fraction, octave_steps) %>%
                     unique(),
                   by = join_by("UUID", `Freq (kHz)` == Frequency)) %>%
         # filter out omitted trials
         group_by(UUID) %>% 
         do(filter(., ! Trial_number %in% .$omit_list)) %>%
         ungroup %>%
         # Apply filters
         filter(detail == "Normal") %>%
         filter(HL_state == "baseline") %>%
         filter(! is.na(octave_steps)) %>%
         # get daily 
         reframe(trial_count = n(),
                 Trial_type = unique(Trial_type),
                 .by = c(rat_ID, rat_name, Genotype, task, detail, 
                         `Freq (kHz)`, octave_fraction, octave_steps)))


## Daily trials ----
# trial_daily_stats = 
  trials_data %>%
  filter(UUID %in% (Discrimination_data$UUID)) %>%
  left_join(Discrimination_data %>%
              select(date, rat_ID, rat_name, Genotype, task, detail, UUID, omit_list) %>%
              unique(),
            by = "UUID") %>%
  # filter out omitted trials
  group_by(UUID) %>% 
  do(filter(., ! Trial_number %in% .$omit_list)) %>%
  ungroup %>%
  # get daily 
  reframe(trials = n(),
          Trial_type = unique(Trial_type),
          Rxn = mean(`Reaction_(s)`, na.rm = TRUE) * 1000,
          delay = mean(`Delay (s)`, na.rm = TRUE),
            .by = c(date, rat_ID, rat_name, Genotype, task, detail, `Freq (kHz)`, Response)) %>%
  mutate(Freq_total_trials = sum(trials),
         .by = c(date, rat_ID, rat_name, Genotype, task, detail, `Freq (kHz)`)) %>%
  mutate(total_trials = sum(trials),
         .by = c(date, rat_ID, rat_name, Genotype, task, detail))


# Descriptive Stats -------------------------------------------------------
cat("getting stats...")
stats_table =
  dataset %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, Genotype, task, detail) %>%
  # Get Averages
  summarise(trial_count = mean(trial_count, na.rm = TRUE),
            hit_percent = mean(hit_percent, na.rm = TRUE),
            FA_percent = mean(FA_percent, na.rm = TRUE),
            .groups = "drop")

names(stats_table)

## Stats ----

### Overall ANOVA ----
discriptive_stats_aov = aov(value ~ Genotype * task * stat,
                            data = stats_table %>%
                              filter(detail == "Normal") %>%
                              filter(task != "Reset") %>%
                              gather(key = "stat", value = "value", c(trial_count:FA_percent)))

Parametric_Check(discriptive_stats_aov)

summary(discriptive_stats_aov)

broom::tidy(TukeyHSD(discriptive_stats_aov)) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value)) %>%
  mutate(Comp1 = str_split_fixed(.$contrast, '-', 2)[,1],
         Comp2 = str_split_fixed(.$contrast, '-', 2)[,2],
         geno1 = str_split_fixed(Comp1, ':', 3)[,1],
         octave_step1 = str_split_fixed(Comp1, ':', 3)[,2] %>% as.numeric(),
         geno2 = str_split_fixed(Comp2, ':', 3)[,1],
         octave_step2 = str_split_fixed(Comp2, ':', 3)[,2] %>% as.numeric(),)

### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "stat", "task" # Main effects
), 
function(x) kruskal.test(reformulate(x, "value"),
                         data = stats_table %>%
                           filter(detail == "Normal") %>%
                           filter(task != "Reset") %>%
                           gather(key = "stat", value = "value", c(trial_count:FA_percent)))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig) %>%
  View

### Post-Hoc Dunn's Test ----
FSA::dunnTest(value ~ interaction(Genotype, stat, task),
              data = stats_table %>%
                filter(detail == "Normal") %>%
                filter(task != "Reset") %>%
                gather(key = "stat", value = "value", c(trial_count:FA_percent)),
              method = "bonf") %>%
  .$res %>%
  as_tibble() %>%
  select(-P.unadj) %>%
  mutate(Sig = gtools::stars.pval(P.adj),
         Comp1 = str_split_fixed(.$Comparison, ' - ', 2)[,1],
         Comp2 = str_split_fixed(.$Comparison, ' - ', 2)[,2],
         geno1 = str_split_fixed(Comp1, '\\.', 3)[,1],
         stat1 = str_split_fixed(Comp1, '\\.', 3)[,2],
         task1 = str_split_fixed(Comp1, '\\.', 3)[,3],
         geno2 = str_split_fixed(Comp2, '\\.', 3)[,1],
         stat2 = str_split_fixed(Comp2, '\\.', 3)[,2],
         task2 = str_split_fixed(Comp1, '\\.', 3)[,3]) %>%
  # only compare within a sex (sib-sib direct comparison)
  filter(stat1 == stat2) %>%
  arrange(stat1) %>%
  # filter(! Sig %in% c(" ", ".")) %>%
  select(-Comp1, -Comp2, -stat2) %>%
  relocate(Comparison, .after = task2) %>%
  View

### trials ----
Parametric_Check(aov(trial_count ~ Genotype * task, 
                     data = stats_table %>%
                       filter(detail == "Normal") %>%
                       filter(task != "Reset")))

summary(aov(trial_count ~ Genotype * task, 
            data = stats_table %>%
              filter(detail == "Normal") %>%
              filter(task != "Reset")))

Parametric_Check(aov(trial_count ~ Genotype, data = filter(stats_table, detail == "Normal")))

trials = t.test(trial_count ~ Genotype, data = filter(stats_table, detail == "Normal"))

### hit % ----
Parametric_Check(aov(hit_percent ~ Genotype * task, 
                     data = stats_table %>%
                       filter(detail == "Normal") %>%
                       filter(task != "Reset")))

#### Non-parametric ----
# Kruskal Testing - Main effects only 
lapply(c("Genotype", "task" # Main effects
), 
function(x) kruskal.test(reformulate(x, "hit_percent"),
                         data = stats_table %>%
                           filter(detail == "Normal") %>%
                           filter(task != "Reset"))) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, p.value, adj.p.value, sig)

Parametric_Check(aov(hit_percent ~ Genotype, data = filter(stats_table, detail == "Normal")))

hit = wilcox.test(hit_percent ~ Genotype, data = filter(stats_table, detail == "Normal"))

### FA % ----
Parametric_Check(aov(FA_percent ~ Genotype * task, 
                     data = stats_table %>%
                       filter(detail == "Normal") %>%
                       filter(task != "Reset")))

summary(aov(FA_percent ~ Genotype * task, 
            data = stats_table %>%
              filter(detail == "Normal") %>%
              filter(task != "Reset")))

Parametric_Check(aov(FA_percent ~ Genotype, data = filter(stats_table, detail == "Normal")))

FA = t.test(FA_percent ~ Genotype, data = filter(stats_table, detail == "Normal"))


## Graph General Stats ----
stats_plot = 
stats_table %>%
  filter(detail == "Normal") %>%
  filter(task != "Reset") %>%
  gather(key = "stat", value = "value", c(trial_count:FA_percent)) %>%
ggplot(aes(x = Genotype, y = value, fill = Genotype, group = Genotype)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       fill = "Genotype") +
  facet_wrap(task ~ stat, scales = "free") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

print(stats_plot)

ggsave(filename = "Fmr1 SD Rxn for tone discrimination.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 150)

# Threshold ---------------------------------------------------------------
# relative to octave steps
cat("calulating thresholds...\n")

Calculate_TH_Octave <- function(df) {
  fit = loess(dprime ~ octave_steps, data = df)
  # plot(fit)
  TH = approx(x = fit$fitted, y = fit$x, xout = TH_cutoff, ties = "ordered")$y
  # print(TH)
  return(TH)
}

octave_TH_table =
  Discrimination_data %>%
  # Sort for ordered (start w/ ID which is unique)
  arrange(rat_ID, rat_name, detail, octave_steps) %>%
  # Prep for Calculate_TH function
  nest(data = c(octave_steps, Frequency, dprime), .by = c(rat_ID, rat_name, Genotype, detail)) %>%
  # mutate not summarise to keep the other columns
  mutate(TH = map_dbl(data, Calculate_TH_Octave)) %>%
  select(-data) %>%
  drop_na()

octave_TH_table_Zoom =
  Discrimination_data %>%
  # Sort for ordered (start w/ ID which is unique)
  arrange(rat_ID, rat_name, detail, octave_steps) %>%
  # Prep for Calculate_TH function
  nest(data = c(octave_steps, Frequency, dprime), .by = c(rat_ID, rat_name, Genotype, detail, Type)) %>%
  # mutate not summarise to keep the other columns
  mutate(TH = map_dbl(data, Calculate_TH_Octave)) %>%
  select(-data) %>%
  drop_na()

TH = 
  octave_TH_table %>%
    reframe(TH = mean(TH, na.rm = TRUE),
            .by = c(Genotype, detail)) %>%
    arrange(detail)

print(TH)

## TH stats -----
octave_TH_table %>%
  filter(detail == "Normal") %>%
  group_by(Genotype) %>%
  do(shapiro.test(.$TH) %>% tidy()) %>%
  mutate(sig = stars.pval(p.value)) %>%
  print

var.test(TH ~ Genotype, 
         data = as.data.frame(octave_TH_table %>% filter(detail == "Normal"))) %>% 
  tidy %>% select(method, p.value) %>% mutate(sig = stars.pval(p.value)) %>% print

t.test(TH ~ Genotype, paired = FALSE, alternative = "two.sided", 
       data = as.data.frame(octave_TH_table %>% filter(detail == "Normal"))) %>% 
  tidy %>% select(method, p.value) %>% mutate(sig = stars.pval(p.value)) %>% print

## Graph TH -----
TH_plot = 
octave_TH_table %>%
  filter(detail == "Normal") %>%
ggplot(aes(x = Genotype, y = TH, fill = Genotype, group = Genotype)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       y = "Threshold (octave step)",
       fill = "Genotype") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

print(TH_plot)

# FA % --------------------------------------------------------------------

discrimination_FA_table =
  Discrimination_data %>%
  filter(detail == "Normal") %>%
  group_by(rat_ID, rat_name, Genotype, octave_steps, Type) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop")

## FA by Freq stats----

FA_freq_data = 
  Discrimination_data %>%
  filter(detail == "Normal") %>%
  # filter(Type == "Broad") %>%
  group_by(rat_ID, rat_name, Genotype, Go_freq, NoGo_freq, detail, octave_steps, Type) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Go_freq = as.factor(Go_freq))

FA_freq_aov = aov(FA_percent_detailed ~ Go_freq * Genotype,
                       data = FA_freq_data)

shapiro.test(FA_freq_aov$residuals)

summary(FA_freq_aov)

kruskal.test(FA_percent_detailed ~ Go_freq,
             data = FA_freq_data)

# Kruskal Testing - Main effects only 
lapply(c("Genotype", "Go_freq", # Main effects
         "interaction(Genotype, Go_freq)"
), 
function(x) kruskal.test(reformulate(x, "FA_percent_detailed"),
                         data = FA_freq_data)) %>% 
  # Convert to table
  do.call(rbind, .) %>% as_tibble() %>% mutate_all(unlist) %>%
  # do a p adjustment and then sig label
  mutate(adj.p.value = p.adjust(p.value, "bonf"),
         sig = gtools::stars.pval(adj.p.value)) %>%
  select(method, parameter, statistic, data.name, adj.p.value, sig)



## Graph of FA -----
FA_plot =
ggplot(data = discrimination_FA_table, 
       aes(x = octave_steps, y = FA_percent_detailed * 100,
           color = Genotype, shape = Type, linetype = Type, group = interaction(Genotype, Type))) +
  geom_hline(yintercept = 50, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  # geom_smooth(se = FALSE) +
  # add the TH line
  # geom_vline(data = TH %>% filter(detail == "Normal"), aes(xintercept = TH, color = Genotype), show.legend = FALSE) +
  # mean for each frequency by genotype
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  # add labels by x-axis
  # geom_text(data = tibble(octave_steps = 12, FA_percent_detailed = 0,
  #                         tone = "No Go", Genotype = "WT"),
  #           aes(label = tone), size = 5, show.legend = FALSE, vjust = 1,
  #           family = "EconSansCndReg") +
  # geom_text(data = tibble(octave_steps = 1, FA_percent_detailed = 0,
  #                         tone = "Go", Genotype = "WT"),
  #           aes(label = tone), size = 4, show.legend = FALSE, vjust = 1,
  #           family = "EconSansCndReg") +
  scale_x_continuous(breaks = c(1, seq(0, 12, by = 2))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       # title = "Discrimination across an octave",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

print(FA_plot)

ggsave(filename = "Fmr1 SD FA for tone discrimination.svg",
       path = save_folder,
       plot = FA_plot,
       width = 10, height = 8, units = "in", dpi = 150)

## FA by Frequency ----
Discrimination_data %>%
  filter(detail == "Normal") %>%
  filter(Type == "Broad") %>%
  group_by(rat_ID, rat_name, Genotype, Go_freq, NoGo_freq, detail, octave_steps, Type) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Go_freq = as.factor(Go_freq),
         NoGo_freq = as.factor(NoGo_freq)) %>%
  ggplot(aes(x = octave_steps, y = FA_percent_detailed * 100,
             color = NoGo_freq, shape = detail, linetype = Genotype, group = interaction(Genotype, NoGo_freq, Type, detail))) +
  geom_hline(yintercept = 50, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  # mean for each frequency by genotype
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2, size = 3) +
  scale_x_continuous(breaks = c(1, seq(0, 12, by = 2))) +
  scale_y_continuous(limits = c(0, 100)) +
  # scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  facet_wrap(~ NoGo_freq) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       linetype = "Genotype",
       shape = "Detail",
       color = "No Go Frequency") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

## Holding FA & hit rate ----
dataset %>%
  filter(detail == "Normal") %>%
  filter(task == "Holding") %>%
  summarise(hit = mean(hit_percent, na.rm = TRUE),
            FA = mean(FA_percent, na.rm = TRUE),
            .by = c(rat_ID, rat_name, Genotype))


# over time ------------------------------------------------------------
discrimination_run_count =
  Discrimination_data %>%
  filter(detail == "Normal") %>%
  filter(Type == "Broad") %>%
  select(date, rat_name, rat_ID, file_name) %>%
  unique %>%
  group_by(rat_ID) %>%
  do(
    arrange(., date) %>%
      rowid_to_column()
  ) %>%
  rename(Measure = rowid)
  
## FA ----
Discrimination_data %>%
  filter(detail == "Normal") %>%
  filter(Type == "Broad") %>%
  left_join(discrimination_run_count,
            by = join_by(date, rat_name, rat_ID, file_name)) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            .by = c(Measure, rat_ID, Genotype, Go_freq, NoGo_freq, detail, octave_steps, Type)) %>%
  mutate(Go_freq = as.factor(Go_freq)) %>%
  filter(Measure < 5) %>%
  mutate(Measure = as.factor(Measure)) %>%
ggplot(aes(x = octave_steps, y = FA_percent_detailed * 100,
           color = Genotype, shape = Measure, linetype = Measure, group = interaction(Genotype, Measure))) +
  geom_hline(yintercept = 50, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  # stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype, size = Measure),
  #              show.legend = TRUE, position = position_dodge(1), vjust = 0.5) +
  scale_x_continuous(breaks = c(1, seq(0, 12, by = 2))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  facet_wrap(~ Genotype) +
  labs(x = "Octave Step",
       y = "False Alarm %",
       # title = "Discrimination across an octave",
       linetype = "Run #", shape = "Run #",
       color = "Genotype") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

ggsave(filename = "Fmr1 SD FA over time for tone discrimination.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 150)

## Hit Rate ----

Discrimination_data %>%
  filter(detail == "Normal") %>%
  filter(Type == "Broad") %>%
  left_join(discrimination_run_count,
            by = join_by(date, rat_name, rat_ID, file_name)) %>%
  summarise(hit_percent = mean(hit_percent, na.rm = TRUE) * 100,
            .by = c(Measure, rat_ID, Genotype, Go_freq, NoGo_freq, detail, Type)) %>%
  mutate(Go_freq = as.factor(Go_freq)) %>%
  filter(Measure < 5) %>%
  mutate(Measure = as.factor(Measure)) %>%
  ggplot(aes(x = Measure, y = hit_percent, fill = Genotype, group = interaction(Measure, Genotype))) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  ylim(75, 100) +
  labs(x = "Run #", y = "Hit %",
       fill = "Genotype") +
  facet_wrap( ~ Genotype, strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Fmr1 SD Hit over time for tone discrimination.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 8, height = 10, units = "in", dpi = 150)

## Rxn Rate ----

Discrimination_data %>%
  filter(detail == "Normal") %>%
  filter(Type == "Broad") %>%
  left_join(discrimination_run_count,
            by = join_by(date, rat_name, rat_ID, file_name)) %>%
  summarise(reaction = mean(reaction, na.rm = TRUE),
            .by = c(Measure, rat_ID, Genotype, Go_freq, NoGo_freq, detail, Type)) %>%
  mutate(Go_freq = as.factor(Go_freq)) %>%
  filter(Measure < 5) %>%
  mutate(Measure = as.factor(Measure)) %>%
  ggplot(aes(x = Measure, y = reaction, fill = Genotype, group = interaction(Measure, Genotype))) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  ylim(100, 500) +
  labs(x = "Run #", y = "Reaction time (Hits only)",
       fill = "Genotype") +
  facet_wrap( ~ Genotype, strip.position = "bottom") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Fmr1 SD Rxn over time for tone discrimination.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 8, height = 10, units = "in", dpi = 150)


# dprime ------------------------------------------------------------------

## graph ====
dprime_plot =
Discrimination_data %>%
  filter(detail == "Normal") %>%
  group_by(rat_ID, rat_name, Genotype, octave_steps, Type) %>%
  summarise(FA_percent_detailed = mean(FA_percent_detailed, na.rm = TRUE),
            dprime_detailed = mean(dprime, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = octave_steps, y = dprime_detailed,
             color = Genotype, shape = Type, linetype = Type, group = interaction(Genotype, Type))) +
  geom_hline(yintercept = 2, color = "forestgreen") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  # mean for genotypes across all frequencies
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "line", linewidth = 1.5, position = position_dodge(.1)) +
  # mean for each frequency by genotype
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               geom = "point", position = position_dodge(.1), stroke = 2) +
  scale_x_continuous(breaks = c(1, seq(0, 12, by = 2))) +
  #scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  labs(x = "Octave Step",
       y = "Sensitivity (d')",
       # title = "Discrimination across an octave",
       fill = "Line", shape = "Line",
       color = "Genotype") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

print(dprime_plot)

ggsave(filename = "Fmr1 SD d' for tone discrimination.svg",
       path = save_folder,
       plot = dprime_plot,
       width = 10, height = 8, units = "in", dpi = 150)


# Reaction time -----------------------------------------------------------

## Averages ----
Reaction_data %>%
  filter(detail == "Normal") %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          .by = c(rat_ID, detail, Genotype)) %>%
  reframe(n = n(),
          reaction_avg = mean(reaction, na.rm = TRUE),
          SD = sd(reaction, na.rm = TRUE),
          SE = FSA::se(reaction, na.rm = TRUE),
          .by = c(Genotype))
  


## Holding Performance ----
dataset %>% 
  filter(task == "Holding") %>%
  filter(detail == "Normal") %>%
  # Select 1st 5 days on holding
    group_by(rat_ID, rat_name, task, detail) %>%
    do(arrange(., desc(date)) %>%
         # Select most recent days for each frequency and task
         tail(n = 5))  %>%
    ungroup %>%
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          n = n(),
          .by = c(rat_ID, rat_name, detail, Genotype)) %>%
  reframe(reaction_avg = mean(reaction, na.rm = TRUE),
          SD = sd(reaction, na.rm = TRUE),
          SE = FSA::se(reaction, na.rm = TRUE),
          .by = c(Genotype))

## Rxn stats -----

### 2 days after training ----
Reaction_2days = 
  dataset %>% 
  filter(task == "Holding") %>%
  filter(detail == "Normal") %>%
  # Select 1st 5 days on holding
  group_by(rat_ID, rat_name, task, detail) %>%
  do(arrange(., desc(date)) %>%
       # Select most recent days for each frequency and task
       head(n = 2))  %>%
  ungroup %>%
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          n = n(),
          .by = c(rat_ID, rat_name, detail, Genotype))

Parametric_Check(aov(reaction ~ Genotype, data = Reaction_2days))

t.test(reaction ~ Genotype, data = Reaction_2days)

### 5 days before reverse ----

Reaction_5days = 
  dataset %>% 
  filter(task == "Holding") %>%
  filter(detail == "Normal") %>%
  # Select 1st 5 days on holding
  group_by(rat_ID, rat_name, task, detail) %>%
  do(arrange(., desc(date)) %>%
       # Select most recent days for each frequency and task
       tail(n = 5))  %>%
  ungroup %>%
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          n = n(),
          .by = c(rat_ID, rat_name, detail, Genotype))

Parametric_Check(aov(reaction ~ Genotype, data = Reaction_5days))

t.test(reaction ~ Genotype, data = Reaction_5days)

### Discrimination days -----
Reaction_discrim = 
  dataset %>% 
  filter(task == "Discrimination") %>%
  filter(detail == "Normal") %>%
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          n = n(),
          .by = c(rat_ID, rat_name, detail, Genotype))

Parametric_Check(aov(reaction ~ Genotype, data = Reaction_discrim))

t.test(reaction ~ Genotype, data = Reaction_discrim)

### Old stats ----

Reaction_data = 
  Discrimination_data %>%
  # remove any duplicates caused by unnesting
  select(rat_ID, reaction, Type, detail, Genotype, UUID) %>% unique %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          .by = c(rat_ID, detail, Genotype, Type)) %>%
  filter(Type != "Error")

Reaction.aov.data = Reaction_data %>%
  filter(detail == "Normal")

Reaction.aov.data$Gaus = LambertW::Gaussianize(Reaction.aov.data$reaction)[, 1]

Rxn.aov = aov(reaction ~ Genotype * Type, data = Reaction.aov.data)

Parametric_Check(Rxn.aov)

summary(Rxn.aov)

broom::tidy(TukeyHSD(Rxn.aov)) %>% 
  mutate(sig = gtools::stars.pval(adj.p.value))

## Rxn Graph ----
Reaction_data %>%
  filter(detail == "Normal") %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          .by = c(rat_ID, detail, Genotype)) %>%
  ggplot(aes(x = Genotype, y = reaction, fill = Genotype, group = Genotype)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       y = "Reaction time (ms)",
       fill = "Genotype") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Fmr1 SD Rxn for tone discrimination.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 150)

## Rxn Graph by file type ----
Reaction_data %>%
  filter(detail == "Normal") %>%
  ggplot(aes(x = Type, y = reaction, fill = Genotype, group = Type)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       y = "Reaction time (ms)",
       fill = "Genotype") +
  facet_wrap(~ Genotype) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Rxn by freq graph ----
Discrimination_data %>%
  filter(detail == "Normal") %>%
  filter(Type == "Broad") %>%
  mutate(Go_freq = as.factor(Go_freq),
         NoGo_freq = as.factor(NoGo_freq)) %>%
  # remove any duplicates caused by unnesting
  select(rat_ID, reaction, Type, detail, Genotype, Go_freq, UUID) %>% unique %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          .by = c(rat_ID, detail, Genotype, Type, Go_freq)) %>%
  filter(Type != "Error") %>%
  ggplot(aes(x = Genotype, y = reaction, fill = Genotype, group = interaction(Genotype))) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       y = "Reaction time (ms)",
       fill = "Genotype") +
  facet_wrap(~ Go_freq) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Holding Rxn ----
dataset %>% 
  filter(task == "Holding") %>%
  filter(detail == "Normal") %>%
  # Select 1st 5 days on holding
  group_by(rat_ID, rat_name, task, detail) %>%
  do(arrange(., desc(date)) %>%
       # Select most recent days for each frequency and task
       tail(n = 5))  %>%
  ungroup %>%
  mutate(reaction = map_dbl(reaction, pluck, "Rxn")*1000) %>%
  reframe(reaction = mean(reaction, na.rm = TRUE),
          n = n(),
          .by = c(rat_ID, rat_name, detail, Genotype)) %>%
  ggplot(aes(x = Genotype, y = reaction, fill = Genotype, group = Genotype)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       y = "Reaction time (ms)",
       fill = "Genotype") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

ggsave(filename = "Fmr1 SD Rxn for holding 5 days before reversal.svg",
       path = "C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Progress Report Figs",
       plot = last_plot(),
       width = 10, height = 8, units = "in", dpi = 150)

# Training stuff ----------------------------------------------------------

## Days in training -----
training_days_data = 
  Training_data %>%
  summarise(days = n(), Genotype = unique(Genotype),
            .by = c(rat_ID, detail))

summarise(training_days_data,
          Days_average = mean(days, na.rm = TRUE),
          SD = sd(days, na.rm = TRUE),
          SE = FSA::se(days, na.rm = TRUE),
          .by = c(Genotype, detail))

t.test(days ~ Genotype,
       data = Training_data %>%
         summarise(days = n(), Genotype = unique(Genotype),
                   .by = c(rat_ID, detail)) %>%
         filter(detail == "Normal"))

training_days_data$Gaus = LambertW::Gaussianize(training_days_data$days)[, 1]

Taining_days.aov = aov(days ~ Genotype * detail, data = training_days_data)

Parametric_Check(Taining_days.aov)

summary(Taining_days.aov)

TukeyHSD(Taining_days.aov) %>% tidy %>%
  select(term, contrast, adj.p.value) %>% 
  mutate(sig = stars.pval(adj.p.value)) %>% print

## Days in each phase/file -----
training_phase_data = 
  Training_data %>%
  summarise(days = n(), Genotype = unique(Genotype), 
            file_name = unique(file_name), date = min(date),
            .by = c(rat_ID, detail, file_name)) %>%
  arrange(rat_ID, date, detail)

### Days Graph ----
training_days_data %>%
  # filter(detail == "Normal") %>%
  ggplot(aes(x = Genotype, y = days, fill = Genotype, group = Genotype)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  stat_summary(fun.data = n_fun, geom = "text", aes(color = Genotype),
               show.legend = FALSE, position = position_dodge(1), vjust = 2, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_fill_manual(values = c("WT" = "darkgrey", "KO" = "red")) +
  labs(x = "",
       y = "Days to criterion",
       fill = "Genotype") +
  facet_wrap(~ detail) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

## Reversal FA & hit -----
Reversal_training = 
  Training_data %>%
  filter(detail == "Reversed") %>%
  # calculate relative date from start of reversal
  group_by(rat_ID) %>%
  do(arrange(., date) %>% mutate(day = row_number()))
  
n_reversal_data =
  Reversal_training %>% ungroup %>% 
  reframe(n = paste("n =", n()), dprime_avg = mean(dprime, na.rm = TRUE),
          .by = c(day, Genotype)) %>%
  mutate(dprime = if_else(Genotype == "WT", dprime_avg + 0.55, dprime_avg - 0.55)) %>%
  group_by(Genotype, n) %>%
  do(arrange(., day) %>% head(., n = 1)) %>% arrange(day) %>% print
  
Reversal_training %>%  
  ggplot(aes(x = day, y = dprime, color = Genotype, group = Genotype)) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE),
               fun.min = function(x) mean(x, na.rm = TRUE) - FSA::se(x, na.rm = TRUE),
               fun.max = function(x) mean(x, na.rm = TRUE) + FSA::se(x, na.rm = TRUE),
               geom = "errorbar", width = 0, position = position_dodge(0.2)) +
  # stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
  #              geom = "line", linewidth = 1.5, position = position_dodge(.2)) +
  geom_smooth(se = FALSE) +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), 
               geom = "point", size = 2, position = position_dodge(.2)) +
  geom_text(data = n_reversal_data, aes(label = n), 
            show.legend = FALSE, size = 3) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_x_continuous(limits = c(0,25), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  ) 

