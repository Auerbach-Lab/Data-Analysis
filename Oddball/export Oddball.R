# Write out for Walker ----------------------------------------------------

AC_Model_data_by_position %>%
  # filter(task == "Base case") %>%
  mutate(detail = factor(detail, 
                         levels = c("Round 1", "Round 2", "Round 3", "Between Treatment", "CNO 3mg/kg")),
         Genotype = str_remove(Genotype, pattern = "Fmr1-LE_")) %>%
  fwrite(glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Walker/Fmr1_AC_CNO_Oddball_Grp1_",
              str_remove_all(Sys.Date(), "-"),".csv"), row.names = FALSE)

table2 =
  AC_Model_data_by_position %>%
  # filter(task == "Base case") %>%
  mutate(detail = factor(detail, 
                         levels = c("Round 1", "Round 2", "Round 3", "Between Treatment", "CNO 3mg/kg")),
         Genotype = str_remove(Genotype, pattern = "Fmr1-LE_")) %>%
  reframe(Rxn = mean(Rxn, na.rm = TRUE),
          .by = c(rat_ID, rat_name, Genotype, Sex, 
                  task, detail, go, Position)) %>%
  group_by(rat_ID, rat_name, Sex, Genotype, task, detail) %>%
  do(mutate(., Rxn_norm = Rxn/filter(., Position == min(Position))$Rxn,
            Rxn_diff = Rxn - filter(., Position == min(Position))$Rxn)) %>%
  ungroup
  
  
  
  
  fwrite(glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Walker/Fmr1_AC_CNO_Oddball_Grp1_",
              str_remove_all(Sys.Date(), "-"),".csv"), row.names = FALSE)
