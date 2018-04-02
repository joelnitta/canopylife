# correlated_evo

# run Pagel's (1994) test of correlated evolution on two binary traits using fitPagel() in phytools

# load packages
library(phytools)
library(mooreaferns)
library(tidyverse)

# load and wrangle data ---------------------------------------------------

# load tree
phy <- mooreaferns::fern_tree

# load un-transformed species traits
traits <- mooreaferns::fern_traits

# make morphotype binary
traits <- mutate (traits, morphotype = case_when (
  morphotype == "cordate" ~ 1,
  morphotype != "cordate" ~ 0))

# run Pagel’s test of correlated evolution --------------------------------

### make function to trim data and run fitPagel() for a trait of interest vs. growth habit
run_fitPagel <- function (selected_trait, traits, phy) {
  
  # trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_trim <- traits[!is.na(traits[selected_trait]), ]
  phy_trim <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% traits_trim$species)])
  traits_trim <- traits_trim[match(traits_trim$species, phy_trim$tip.label), ]
  
  # extract named vector of trait values
  trait_vec <- traits_trim[[selected_trait]]
  names(trait_vec) <- traits_trim$species
  
  # extract named vector of growth habit
  habit_vec <- traits_trim[["habit"]]
  names(habit_vec) <- traits_trim$species
  
  # run fitPagel() on selected trait vs. growth habit
  model <- fitPagel(tree = phy_trim, x = trait_vec, y = habit_vec)
  
  # get model results
  list(trait = selected_trait,
       logL_indep = model$independent.logL, 
       logL_dep = model$dependent.logL, 
       likelihood_ratio = model$lik.ratio, 
       pval = model$P)

}

# set up input for mapping run_fitPagel() to list of selected traits
map_input <- list(
  as.list(c("gemmae", "glands", "hairs", "morphotype")),
  list(traits),
  list(phy)
)

# map run_fitPagel() to selected traits
pagel.results <- pmap_dfr(map_input, run_fitPagel)

# rename rows for xtable
pagel.results <- as.data.frame(pagel.results)
rownames(pagel.results) <- pagel.results$trait
pagel.results$trait <- NULL

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
