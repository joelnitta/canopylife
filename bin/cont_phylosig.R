# cont_phylosig.R

# Analyze phylogenetic signal in continuous traits
# using Blomberg's K and Pagel's lambda

# load packages
library(phytools) # phylosig
library(mooreaferns)
library(tidyverse)

# set working directory
setwd(here::here())

# load data ---------------------------------------------------------------

# tree
phy <- mooreaferns::fern_tree

# un-transformed species traits
traits <- mooreaferns::fern_traits

# Analyze Blomberg's K and Pagel's Lambda ---------------------------------

### make function to trim data and run phylosig() for a trait of interest
run_phylosig <- function (selected_trait, traits, phy) {
  
  # trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_trim <- traits[!is.na(traits[selected_trait]), ]
  phy_trim <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% traits_trim$species)])
  traits_trim <- traits_trim[match(traits_trim$species, phy_trim$tip.label), ]

  # extract named vector of trait values for phylosig()
  trait_vec <- traits_trim[[selected_trait]]
  names(trait_vec) <- traits_trim$species
  
  # run phylosig() on selected trait
  # using Blomberg's K
  k_model <- phylosig(phy_trim, trait_vec, method = "K", test = TRUE)
  # and Pagel's lambda
  lambda_model <- phylosig(phy_trim, trait_vec, method = "lambda", test = TRUE)
  
  # get model results
  list(trait = selected_trait,
       kval = k_model$K,
       k.pval = k_model$P,
       lambda = lambda_model$lambda,
       lambda.pval = lambda_model$P)

}

# set up input for mapping run_phylosig() to list of selected traits
map_input <- list(
  as.list(c("stipe", "length", "width", "dissection", "pinna", "sla", "rhizome")),
  list(traits),
  list(phy)
)

# map run_phylosig() to selected traits
cont_phylosig.results <- pmap_dfr(map_input, run_phylosig)

# rename rows for xtable
cont_phylosig.results <- as.data.frame(cont_phylosig.results)
rownames(cont_phylosig.results) <- cont_phylosig.results$trait
cont_phylosig.results$trait <- NULL

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
