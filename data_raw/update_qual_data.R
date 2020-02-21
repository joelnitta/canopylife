library(assertr)
library(tidyverse)

# Read in list of species used in this study
species_list <- read_csv(
  "data/nitta_2017/species.csv", 
  col_types = "cccccnn") %>%
  janitor::clean_names() %>%
  filter(include == 1, tahiti_only == 0) %>%
  pull(genus_sp)

# Read in previous version of qualitative trait data
morph_qual <- read_csv(
  "data_raw/morph_qual_traits_2020-02-11.csv", 
  na = c("?", "n/a", "NA", "N/A", ""),
  col_types = "cnnnncncc")

# Make a tibble of Prosaptia subnuda data to add
prosaptia_subnuda_data <- tibble(
  species = "Prosaptia_subnuda",
  dissection = "pinnatifid",
  growth_habit = "epiphytic",
  morphotype = "strap",
  source = "4,5"
)

# Update data
morph_qual_updated <- 
morph_qual %>%
  # Filter to only species in species list
  filter(species %in% species_list) %>%
  # Use character instead of numeric for dissection
  mutate(
    dissection = case_when(
      dissection == 1 ~ "simple",
      dissection == 2 ~ "pinnatifid",
      dissection == 3 ~ "one-pinnate",
      dissection == 4 ~ "one-pinnate-pinnatifid",
      dissection == 5 ~ "two-pinnate",
      dissection == 6 ~ "two-pinnate-pinnatifid",
      dissection == 7 ~ "three-pinnate",
      dissection == 8 ~ "three-pinnate-pinnatifid",
      dissection == 9 ~ "more_divided",
      dissection == 10 ~ "bipinnatifid_or_tripinnatifid",
      TRUE ~ NA_character_
    )
  ) %>%
  # Change dissection for some species
  mutate(
    dissection = case_when(
      species == "Asplenium_affine" ~ "two-pinnate-pinnatifid", # from once-pinnate
      species == "Tmesipteris_gracilis" ~ "pinnatifid", # from simple
      species == "Humata_anderssonii" ~ "two-pinnate-pinnatifid", # from once-pinnate-pinnatifid
      species == "Crepidomanes_kurzii" ~ "pinnatifid", # from once-pinnate
      species == "Diplopterygium_longissimum" ~ "two-pinnate", # missing in last version
      species == "Hymenophyllum_braithwaitei" ~ "bipinnatifid_or_tripinnatifid", # missing in last version
      species == "Hypolepis_sp1" ~ "three-pinnate-pinnatifid", # missing in last version
      species == "Lindsaea_propinqua" ~ "two-pinnate", # missing in last version
      TRUE ~ dissection
    )
  ) %>%
  # Add column for growth habit (not only binary)
  assert(not_na, habit_binary) %>%
  assert(in_set(0, 1), habit_binary) %>%
  mutate(
    growth_habit = case_when(
      habit_binary == 0 ~ "terrestrial",
      habit_binary == 1 ~ "epiphytic")
  ) %>%
  # Update growth habit for some species
  mutate(
    growth_habit = case_when(
      species == "Adiantum_hispidulum" ~ "epipetric",
      species == "Adiantum_raddianum" ~ "epipetric",
      species == "Arthropteris_palisotii" ~ "hemiepiphytic",
      species == "Cheilanthes_nudiuscula" ~ "terrestrial",
      species == "Crepidomanes_kurzii" ~ "epipetric",
      species == "Lomagramma_tahitensis" ~ "hemiepiphytic",
      species == "Lomariopsis_brackenridgei" ~ "hemiepiphytic",
      species == "Lygodium_reticulatum" ~ "climbing",
      species == "Teratophyllum_wilkesianum" ~ "climbing",
      TRUE ~ growth_habit
    )
  ) %>%
  # Change gametophyte morphotype for some species (all to strap from ribbon)
  mutate(
    morphotype = case_when(
      species == "Calymmodon_orientalis" ~ "strap",
      species == "Grammitis_marginelloides" ~ "strap",
      species == "Humata_anderssonii" ~ "strap",
      species == "Humata_pectinata" ~ "strap",
      TRUE ~ morphotype
    )
  ) %>%
  # Add formerly missing P. subnuda data
  bind_rows(prosaptia_subnuda_data) %>%
  # Final checks
  assert(not_na, growth_habit) %>%
  assert(in_set("epiphytic", "epipetric", "terrestrial", "hemiepiphytic", "climbing"), growth_habit) %>%
  # Reorder columns
  select(species, growth_habit, dissection, glands, hairs, morphotype, gemmae, source, notes)

write_csv(morph_qual_updated, "data/morph_qual_traits.csv")
