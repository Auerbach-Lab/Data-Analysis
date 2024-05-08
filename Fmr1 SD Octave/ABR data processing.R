# Load Necessary Datasets -------------------------------------------------

# get list of CSVs to append together
files = list.files(ABR_data_folder, pattern = "\\d+.csv$", recursive = TRUE, full.names = TRUE)

load_csv_file <- function(file_name) {
  rat_ID = str_remove(file_name, pattern = ABR_data_folder) %>% str_remove(".csv") %>% as.numeric()
  
  csv_table = fread(file_name)
  
  output = as_tibble(csv_table) %>%
    # drop Wave not picked
    filter(Wave_amp != 0 & Wave_lat != 0) %>%
    mutate(rat_ID = rat_ID)
  
  return(output)
}

# load each CSV file generated from the abr_processor script (credit: Walker Gauthier)
ABR_data = lapply(files, load_csv_file) %>% bind_rows()

rat_genotypes = fread(glue("{ABR_data_folder}/Genotypes.csv"))

ABR_data = left_join(ABR_data, rat_genotypes, by = join_by(rat_ID))



# Write our data for Ben --------------------------------------------------

fwrite(ABR_data, glue("C:/Users/Noelle/Box/Behavior Lab/Shared/Ben/Fmr1_SD_ABR_data", str_remove_all(Sys.Date(), "-"),".csv"), row.names = FALSE)
