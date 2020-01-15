# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

make_ref_list(
  rmd_file = "ms/manuscript.Rmd", 
  raw_bib = "ms/references_raw.bib",
  final_bib = "ms/references.bib",
  strip_fields = "abstract",
  exclude = c("Laliberte2011", "Orme2013", "Thiers2019"))

make_ref_list(
  rmd_file = "ms/SI.Rmd", 
  raw_bib = "ms/references_raw.bib", 
  final_bib = "ms/references_si.bib",
  exclude = "FloraNSW2019"
)

# Make some manual fixes to authors in SI bibliography
# (these are institutions, so need double brackets to
# avoid latex thinking they have first and last names)
read_lines("ms/references_si.bib") %>%
  str_replace(
    "National Museum of Nature and Science",
    "\\{National Museum of Nature and Science\\}") %>%
  write_lines("ms/references_si.bib")

# Make some manual fixes to main bibliography
read_lines("ms/references.bib") %>%
  str_replace(
    "Pteridophyte Phylogeny Group I",
    "\\{Pteridophyte Phylogeny Group I\\}") %>%
  str_replace(
    "Adaptive radiation versus 'radiation' and 'explosive diversification'",
    "Adaptive radiation versus `radiation' and `explosive diversification'") %>%
  str_replace(
    "'Sun leaves' and 'shade leaves'",
    "`Sun leaves' and `shade leaves'") %>%
  write_lines("ms/references.bib")
