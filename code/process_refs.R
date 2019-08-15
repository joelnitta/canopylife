# Clean up latex bibliography exported from Mendeley
#
# .bib files in data_raw were exported from Mendeley (select
# all references -> File:export -> save as .bib).
#
# These cannot be used-as is when compiling Rmd: many are formatted
# incorrectly, and I don't need all the references in my library,
# just those cited in the text. 
# 
# This script selects on the cited ones,
# cleans them, and saves the clean bib file.
#
# Since other people don't need the full raw bibliography, 
# this script is run as needed and the results (the clean bib)
# checked-in to the repo.
# The script is left here though for reproducibility.

library(jntools)
library(tidyverse)
library(bib2df)

setwd(here::here())

#' Produce a filtered, cleaned .bib file
#' 
#' Takes the raw .bib file (e.g., one exported from
#' Mendeley with thousands of entries), filters it
#' to only the citations in the Rmd file, and cleans
#' the references so pandoc won't crash.
#'
#' @param rmd_file Path to rmd file
#' @param raw_bib Path to raw bib file
#' @param final_bib Path to write cleaned, filtered bib file
#'
#' @return Nothing
#' 
make_ref_list <- function(rmd_file, raw_bib, final_bib) {

# Process manuscript Rmd and pull out citations
# (words that begin with '@')
citations <-
  read_lines(rmd_file) %>%
  str_split(" |;") %>%
  unlist %>%
  magrittr::extract(., str_detect(., "@")) %>%
  str_remove_all("\\[|\\]|\\)|\\(|\\.$|,") %>%
  unique %>%
  sort %>%
  str_remove_all("@")

# Read in entire raw bibliography exported from Mendeley as tibble
bib_df <- bib2df(raw_bib)

# Select only needed columns and subset to citations in Rmd
bib_df_selected <-
  bib_df %>%
  select(CATEGORY, BIBTEXKEY, AUTHOR, 
         BOOKTITLE, CHAPTER, EDITOR, EDITION,
         JOURNAL, PAGES, PUBLISHER, TITLE, VOLUME, YEAR) %>%
  filter(BIBTEXKEY %in% citations)

# Fix formatting of "Jr."
bib_df_selected <-
  bib_df_selected %>%
  mutate(AUTHOR = map(AUTHOR, ~str_replace_all(., "^Watkins Jr.\\}", "\\{Watkins Jr.\\}")))

# Fix missing bracket for italics at start of title
bib_df_selected <-
  bib_df_selected %>%
  mutate(TITLE = map(TITLE, ~str_replace_all(., "^textless\\}i", "\\{\\\\textless\\}i")))

# Write this out to temporary file so it can be cleaned
# with jntools::clean_bib(), since that works by
# reading in an external bib flie.
temp_file <- tempfile()
df2bib(bib_df_selected, file = temp_file)

# Do final cleaning with jntools::clean_bib
jntools::clean_bib(temp_file) %>% 
  write_lines(final_bib)

}

make_ref_list("ms/manuscript.Rmd", "ms/references_raw.bib", "ms/references.bib")

make_ref_list("ms/SI.Rmd", "ms/references_raw.bib", "ms/si_references.bib")
