# binary_phylosig.R

# Analyze phylogenetic signal in binary (qualitative) traits
# using Fritz and Purvis' D

# load packages
library(caper) # phylo.d()
library(mooreaferns) # data
library(tidyverse) # bind_rows(), map(). load last to avoid unwanted conflicts.

# set working directory
setwd(here::here())

# load and wrangle data ---------------------------------------------------

# load tree
phy <- mooreaferns::fern_tree

# load un-transformed species traits
traits <- mooreaferns::fern_traits

# make morphotype into a binary numeric category: 0 is noncordate, 1 is cordate
traits <- mutate (traits, morphotype = case_when (
  morphotype == "cordate" ~ 1,
  morphotype != "cordate" ~ 0 ))

# make growth habit into a binary category: 0 is not epiphytic (ie, terrestrial), 1 is epiphytic
traits <- mutate (traits, habit = case_when (
  habit == "terrestrial" ~ 0,
  habit == "epiphytic" ~ 1 ))

# set up data lists for phylo.d -------------------------------------------

# make vector of traits to test
binary_traits <- c("habit", "gemmae", "glands", "hairs", "morphotype")

# make list of dataframes, each only including species and the trait to test
traits_list <- map(binary_traits, function (x) {
  traits[,c("species", x)]
})

names(traits_list) <- binary_traits

### make caper tree/trait datasets silently dropping missing data
# some of the traits include NAs.
# for phylo.d, need to have matching tree and trait data with all NAs removed
# don't do comparative.data on all the traits together, or species missing data for ANY trait will get dropped

comp_data_list <- map(traits_list, function (x) {
   comparative.data(phy=phy, data=x, names.col="species", na.omit=TRUE)
  }
)

# run Fitz and Purvis' D  -------------------------------------------------

# helper function to run phylo.d and store observed D and associated P values for each trait
# use eval(parse(text())) to feed variable into phylo.d
# see https://stackoverflow.com/questions/12516260/defining-dependent-and-independent-variables-dynamically-in-the-ezanova-function

run_phylo_d <- function(binary_var) {
  phylo.d.out <- eval(parse
       (text=paste0('phylo.d(data=comp_data_list[["',binary_var,'"]],
                      binvar=',binary_var,')')
  ))
  list(trait = phylo.d.out$binvar, 
       num_present = phylo.d.out$StatesTable[1], 
       num_absent = phylo.d.out$StatesTable[2], 
       D = phylo.d.out$DEstimate, 
       prob_random = phylo.d.out$Pval1, 
       prob_brownian = phylo.d.out$Pval0)
}

# apply function to list of traits
binary_phylosig.results <- map_df(binary_traits, run_phylo_d)

# format results table ----------------------------------------------------

# need to have results as data frame with named rows for xtable
binary_phylosig.results <- as.data.frame(binary_phylosig.results)
rownames(binary_phylosig.results) <- binary_phylosig.results$trait
binary_phylosig.results$trait <- NULL

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
