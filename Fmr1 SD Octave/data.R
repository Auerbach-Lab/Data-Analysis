# Load Necessary Datasets -------------------------------------------------
load(glue("{projects_folder}/run_archive.Rdata"), .GlobalEnv)
rat_archive = fread(glue("{projects_folder}/rat_archive.csv"), 
                    select = c("Rat_ID", "DOB", "Sex", "Genotype", "HL_date"))

# # Individual Trial Data
# load(paste0(projects_folder, "TTS_archive.Rdata"), .GlobalEnv)


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
  mutate(octave_fraction = log(as.numeric(str_extract(file_name, pattern = "[:digit:]+?(?=-.+?kHz)"))/Frequency)/log(2),
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
       y = "Threshold (octave step)",
       fill = "Genotype") +
  facet_wrap(task ~ stat, scales = "free") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

print(stats_plot)

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

# Reaction time -----------------------------------------------------------

Reaction_data = 
  Discrimination_data %>%
    # remove any duplicates caused by unnesting
    select(rat_ID, reaction, Type, detail, Genotype, UUID) %>% unique %>%
    reframe(reaction = mean(reaction, na.rm = TRUE),
            .by = c(rat_ID, detail, Genotype, Type)) %>%
    filter(Type != "Error")


## Rxn stats -----
Reaction.aov.data = Reaction_data %>%
  filter(detail == "Normal")

Reaction.aov.data$Gaus = LambertW::Gaussianize(Reaction.aov.data$reaction)[, 1]

Rxn.aov = aov(reaction ~ Genotype * Type, data = Reaction.aov.data)

Parametric_Check <- function(AOV.data) {
  is_parametric = shapiro.test(AOV.data$residuals)$p.value > 0.05
  
  if (is_parametric == TRUE) {writeLines("Normal data proced with ANOVA")} else
  {writeLines(paste("Non-parametric data so use Kruskal followed by Dunn testing. \nShapiro Test: p =", shapiro.test(AOV.data$residuals)$p.value %>% round(digits = 3)))}
  
}

Parametric_Check(Rxn.aov)

summary(Rxn.aov)

## Rxn Graph ----
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
       y = "Threshold (octave step)",
       fill = "Genotype") +
  facet_wrap(~ Genotype) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255))
  )

# Training stuff ----------------------------------------------------------

## Days in training -----
training_days_data = 
  Training_data %>%
  summarise(days = n(), Genotype = unique(Genotype),
            .by = c(rat_ID, detail))

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
       y = "Threshold (octave step)",
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
  mutate(dprime = if_else(Genotype == "WT", dprime_avg + 0.55, dprime_avg -0.55)) %>%
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

