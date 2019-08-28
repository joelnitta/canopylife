# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

make_ref_list(
  rmd_file = "ms/manuscript.Rmd", 
  raw_bib = "ms/references_raw.bib",
  final_bib = "ms/references.bib")

make_ref_list(
  rmd_file = "ms/SI.Rmd", 
  raw_bib = "ms/references_raw.bib", 
  final_bib = "ms/si_references.bib"
)

# Make some manual fixes to authors in SI bibliography
# (these are institutions, so need double brackets to
# avoid latex thinking they have first and last names)
read_lines("ms/si_references.bib") %>%
  str_replace(
    "National Museum of Nature and Science",
    "\\{National Museum of Nature and Science\\}") %>%
  str_replace(
    "Royal Botanic Gardens and Domain Trust",
    "\\{Royal Botanic Gardens and Domain Trust\\}") %>%
  write_lines("ms/si_references.bib")
