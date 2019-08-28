# Make clean bib files for each Rmd file that only includes cited references

library(jntools)

make_ref_list("ms/manuscript.Rmd", "ms/references_raw.bib", "ms/references.bib")

make_ref_list("ms/SI.Rmd", "ms/references_raw.bib", "ms/si_references.bib")
