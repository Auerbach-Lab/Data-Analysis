# Get core data -----------------------------------------------------------
source("main.R")

# Calculate Overall TH ----------------------------------------------------
cat("Caluclating: thresholds...")

# Threshold calculation calculation based on TH_cutoff intercept of fit curve
# LOESS: Local Regression is a non-parametric approach that fits multiple regressions
# see http://r-statistics.co/Loess-Regression-With-R.html
Calculate_TH <- function(df) {
  rat_name = unique(df$rat_name)
  Freq = unique(df$Freq)
  Dur = unique(df$Dur)
  step_size = if(is.null(df$step_size)) "All data" else glue("{unique(df$step_size)}dB step size")
  
  ## Default model ----------------------------------------------------------
  # fit = loess(dprime ~ dB, data = df)
  # TH = approx(x = fit$fitted, y = fit$x, xout = TH_cutoff, ties = "ordered")$y
  
  # # Get TH point (it needs inverting for some reason)
  # TH_point = approx(x = fit$fitted, y = fit$x, xout = TH_cutoff, ties = "ordered")
  # TH_x = TH_point$y
  # TH_y = TH_point$x
  
  # # Actually plot
  # plot(fit,
  #      main = glue("{rat_name} @ {Freq}kHz & {Dur}ms {step_size} TH: {round(TH, digits = 1)}"),
  #      sub = glue("TH: x = {TH_x}"))
  # # added fitted line
  # lines(fit$x, predict(fit), col = "green")
  # # add TH point to plot
  # points(TH_x, TH_y, col="red")
  # # add line to plot
  # abline(h=1.5, col="blue")
  # # Save
  # dev.copy(png,
  #          glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Rollout Issues/dprime plots/{rat_name} {Freq}kHz {Dur}ms {step_size}.png"))
  # dev.off()
  
  ## DRDA Model ---------------------------------------------------------------
  # Dose-dependant curve which suggests 4-parameter logistic function is best fit
  fit_drda = drda(dprime ~ dB, data = df)
  # pull the 1st number which is the estimate, the other 2 are the 95% CI range
  TH_drda = effective_dose(fit_drda, y = TH_cutoff, type = "absolute")[,1]
  
  # # Plot DRDA
  # plot(fit_drda,
  #      main = glue("{rat_name} @ {Freq}kHz & {Dur}ms {step_size}, TH: {round(TH_drda, digits = 1)}"))
  # # Save
  # dev.copy(png,
  #          glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Noelle/Rollout Issues/dprime plots/{rat_name} {Freq}kHz {Dur}ms {step_size} ddra.png"))
  # dev.off()
  
  
  # print(glue("{rat_name} @ {Freq}kHz & {Dur}ms: {round(TH, digits = 1)} or DRDA {round(TH_drda, digits = 1)}"))
  return(TH_drda)
}

TH_table = core_data %>%
  # Omit Training & Reset days
  filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # Get dprimes
  unnest(dprime) %>%
  # Sort for ordered
  arrange(rat_ID, rat_name, Freq, Dur, dB) %>%
  #Prep for Calculate_TH function
  nest(data = c(rat_name, Freq, Dur, dB, dprime), .by = c(rat_ID, rat_name, sex, genotype, detail, Freq, Dur)) %>% 
  mutate(TH = map_dbl(data, Calculate_TH)) %>%
  select(-data) 


# Get Reaction times by rat -----------------------------------------------
cat("reaction times...")

Rxn_table = core_data %>%
  # Omit Training & Reset days
  dplyr::filter(! task %in% c("Training", "Reset")) %>%
  # Omit days with > 45% FA, i.e. guessing
  filter(FA_percent < FA_cutoff) %>%
  # Get Reaction times:
  unnest(reaction) %>%
  # Use rat_ID because its sure to be unique
  group_by(rat_ID, rat_name, sex, genotype, detail, `Freq (kHz)`, `Dur (ms)`, `Inten (dB)`) %>%
  # Get Averages
  transmute(Rxn = mean(Rxn, na.rm = TRUE) * 1000) %>% 
  unique()

# Limit to Over TH
TH_filter <- function(df) {
  ID = unique(df$ID)
  Dur = unique(df$Dur)
  kHz = unique(df$Freq)
  # kHz = if_else(kHz == "0", "BBN", paste0(kHz,"kHz")) # %>% print
  
  cuttoff = TH_table %>% # may have to use UQ to force the evaluation of the variable
    filter(Dur == Dur & rat_ID == ID & Freq == kHz) %>% 
    .$TH
  
  cuttoff = ifelse(identical(cuttoff, numeric(0)), -99, cuttoff)
  
  r = df %>%
    filter(Inten >= UQ(cuttoff))
  
  return(r)
}

Rxn_table_over_TH = Rxn_table %>% 
  # Prep for TH_filter function
  ungroup() %>%
  mutate(ID = rat_ID, Dur = `Dur (ms)`, Freq = `Freq (kHz)`, Inten = `Inten (dB)`) %>%
  nest(data = c(ID, Freq, Dur, Inten, Rxn), .by = c(rat_ID, rat_name, sex, genotype, detail, `Freq (kHz)`, `Dur (ms)`)) %>%
  # Apply TH_filter
  mutate(data = map(data, TH_filter)) %>%
  unnest(data) %>%
  select(- any_of(c("Dur (ms)", "Freq (kHz)", "ID"))) %>%
  rename(Intensity = Inten, Duration = Dur, Frequency = Freq) %>%
  # Use rat_ID because its sure to be unquie
  group_by(rat_ID, rat_name) %>%
  # Get Averages
  mutate(Rxn = mean(Rxn, na.rm = TRUE)) %>% 
  unique()


# Time to learning --------------------------------------------------------
# cat("learning times...")
# 
# Learning_streak = core_data %>%
#   dplyr::filter(task %in% c("Training")) %>%
#   group_by(rat_ID) %>% # pregroup this so that numbering restarts at 1 for each rat
#   mutate(groupid = data.table::rleid(phase, task, detail)) %>%
#   group_by(rat_ID, sex, genotype, groupid, phase, task, detail) %>%
#   summarise(running_streak = n(), .groups = "drop")


# Renaming -----------------------------------------------------
cat("updating column names...")

TH_table = rename(TH_table, Frequency = Freq, Duration = Dur)

Rxn_table = rename(Rxn_table, Frequency = `Freq (kHz)`, Duration = `Dur (ms)`, Intensity = `Inten (dB)`)

cat("done.\n")

