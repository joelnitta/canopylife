# This script looks up scientific names (specis names with taxonomic author) using
# a live search of taxonomic databases via the Global Name Resolver 
# (https://resolver.globalnames.org/). Since results have the possibility
# of changing, it is run once and the results
# checked-in to the repo at the time of manuscript submission.
# The script is left here though for reproducibility.

library(taxize)
library(janitor)
library(tidyverse)

# Read in list of accepted species names used in this study,
# extract the genus and specific epithet
species_list <- read_csv("data/nitta_2017/species.csv") %>%
  clean_names() %>%
  filter(include == 1, tahiti_only == 0) %>%
  pull(genus_sp)

# Reformat the genus_species list of names into a dataframe
# with the taxon_code used to refer to that taxon in scripts,
# and other names at various taxonomic ranks.
taxa_names <-
  tibble(taxon_code = species_list) %>%
  separate(taxon_code, c("genus", "specific_epithet"), by = "_", remove = FALSE) %>%
  mutate(informal_variety = str_match(specific_epithet, "[0-9]") %>% map_chr(1)) %>%
  mutate(specific_epithet = str_remove_all(specific_epithet, "[0-9]")) %>%
  mutate(species = paste(genus, specific_epithet)) %>%
  mutate(informal_variety = case_when(
    specific_epithet == "sp" ~ NA_character_,
    TRUE ~ informal_variety
  )) %>%
  mutate(infaspecific_epithet = case_when(
    species == "Lindsaea repens" ~ "marquesensis",
    species == "Davallia solida" ~ "solida",
    species == "Ctenitis sciaphila" ~ "sciaphila",
    species == "Dicranopteris linearis" ~ "linearis",
    species == "Deparia petersenii" ~ "congrua",
    TRUE ~ NA_character_
  )) %>%
  mutate(taxon = jntools::paste3(genus, specific_epithet, infaspecific_epithet))

# Use the Global Name Resolver (https://resolver.globalnames.org/) to match the names
# to either IPNI or TROPICOS databases.
taxa_names_with_authors <- gnr_resolve(
  names = taxa_names$taxon, 
  best_match_only = TRUE, 
  # for data_source_ids, 167 = ipni, 165 = tropicos
  data_source_ids = c(167, 165)) %>%
  select(species = user_supplied_name, scientific_name = matched_name, name_source = data_source_title) %>%
  right_join(taxa_names) %>% 
  select(taxon_code, genus, specific_epithet, infaspecific_epithet, informal_variety, scientific_name, name_source) %>%
  # Manually add names from NCBI database or Murdock and Smith 2003
  # for those missing from IPNI or TROPICOS
  mutate(scientific_name = case_when(
    taxon_code == "Diplazium_grantii" ~ "Diplazium grantii (Copel.) C.Chr.",
    taxon_code == "Dryopteris_macrolepidota" ~ "Dryopteris macrolepidota Copel.",
    taxon_code == "Psilotum_nudum" ~ "Psilotum nudum (L.) P.Beauv.",
    taxon_code == "Archigrammitis_tahitensis" ~ "Archigrammitis tahitensis (C.Chr.) Parris",
    taxon_code == "Hymenophyllum_polyanthos" ~ "Hymenophyllum polyanthos (Sw.) Sw.",
    taxon_code == "Ctenitis_sciaphila" ~ "Ctenitis sciaphila (Maxon) Ching var. sciaphila",
    taxon_code == "Davallia_solida" ~ "Davallia solida (G. Forst.) Sw. var. solida",
    taxon_code == "Deparia_petersenii" ~ "Deparia petersenii (Kunze) M. Kato subsp. congrua (Brack.) M. Kato",
    taxon_code == "Dicranopteris_linearis" ~ "Dicranopteris linearis (Burm. f.) Underw. var. linearis",
    taxon_code == "Lindsaea_repens" ~ "Lindsaea repens (Bory) Thwaites var. marquesensis E. D. Br.",
    TRUE ~ scientific_name
  )) %>%
  mutate(name_source = case_when(
    taxon_code == "Diplazium_grantii" ~ "NCBI",
    taxon_code == "Dryopteris_macrolepidota" ~ "NCBI",
    taxon_code == "Psilotum_nudum" ~ "NCBI",
    taxon_code == "Archigrammitis_tahitensis" ~ "NCBI",
    taxon_code == "Hymenophyllum_polyanthos" ~ "NCBI",
    taxon_code == "Ctenitis_sciaphila" ~ "Murdock and Smith 2003",
    taxon_code == "Davallia_solida" ~ "Murdock and Smith 2003",
    taxon_code == "Deparia_petersenii" ~ "Murdock and Smith 2003",
    taxon_code == "Dicranopteris_linearis" ~ "Murdock and Smith 2003",
    taxon_code == "Lindsaea_repens" ~ "Murdock and Smith 2003",
    TRUE ~ name_source
  )) %>%
  mutate(name_source = case_when(
    name_source == "Tropicos - Missouri Botanical Garden" ~ "TROPICOS",
    name_source == "The International Plant Names Index" ~ "IPNI",
    TRUE ~ name_source
  )) %>%
  rename_all(~str_replace(., "_", " "))
  
write_csv(taxa_names_with_authors, here::here("ms/Table_S1.csv"))