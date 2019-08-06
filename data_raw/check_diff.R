# Check for differences in morph_qual_traits_2016-09-03_mod.csv from Moorea_ferns_morph_2016-09-03.xlsx
# I know that I copied morph_qual_traits_2016-09-03_mod.csv from the Moorea_ferns_morph_2016-09-03.xlsx
# morphology worksheet, but didn't document further changes from there at the time.
# This script compares them to find where I changed things.

library(tidyverse)
library(readxl)

# Moorea_ferns_morph_2016-09-03.xlsx is read-only and
# copied from /Volumes/Transcend/Projects/Moorea CommPhy/Morphology
traits_raw <- read_excel("data_raw/Moorea_ferns_morph_2016-09-03.xlsx", sheet = "morphology", skip = 11)
traits_new <- read_csv("data_raw/morph_qual_traits_2016-09-03_mod.csv", na = c("#DIV/0!", "#N/A", ""))

# The order had been rearranged in excel, so sort both by OTU
traits_raw <- traits_raw %>% arrange(OTU)
traits_new <- traits_new %>% arrange(OTU)

# Not the same yet, due to differences in the `Key Group` column
all.equal(traits_raw, traits_new)

# Not using that column any more anyways, so drop it
traits_raw <- traits_raw %>% select(-`Key Group`)
traits_new <- traits_new %>% select(-`Key Group`)

# Differences due to included taxa (rows)
all.equal(traits_raw, traits_new)

# Outgroups have been dropped from the new dataset
setdiff(traits_new$OTU, traits_raw$OTU)
setdiff(traits_raw$OTU, traits_new$OTU)

# Drop these and compare again
traits_raw <- traits_raw %>% filter(traits_raw$OTU %in% traits_new$OTU)

# Most remaining differences are due to coding of NA values
all.equal(as.data.frame(traits_raw), as.data.frame(traits_new), check.attributes = FALSE)

left_join(
select(traits_raw, OTU, dissection_raw = dissection),
select(traits_new, OTU, dissection_new = dissection)
) %>%
  mutate_at(vars(dissection_raw, dissection_new), ~replace_na(., "NA")) %>%
  filter(dissection_raw != dissection_new)

left_join(
  select(traits_raw, OTU, hairs_raw = hairs),
  select(traits_new, OTU, hairs_new = hairs)
) %>%
  mutate_at(vars(hairs_raw, hairs_new), ~replace_na(., "NA")) %>%
  filter(hairs_raw != hairs_new)

left_join(
  select(traits_raw, OTU, gemmae_raw = gemmae),
  select(traits_new, OTU, gemmae_new = gemmae)
) %>%
  mutate_at(vars(gemmae_raw, gemmae_new), ~replace_na(., "NA")) %>%
  filter(gemmae_raw != gemmae_new)
