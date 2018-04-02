# binaryPGLMM.R
# Construct binary phylogenetic generalized linear mixed models (PGLMMs) 
# for binary (gametophyte) traits related to epiphytic growth

# load packages
library(ape) # binaryPGLMM()
library(mooreaferns)
library(tidyverse)

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

# binaryPGLMM -------------------------------------------------------------

### make function to trim data and run binaryPGLMM
# binaryPGLMM can't take NA values, so trim trait dataframes and 
# phylogenies so they only include species with no NA values for each of 
# traits to test
run_bPGLMM <- function (selected_trait, traits, phy) {
  
  # trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_trim <- traits[!is.na(traits[selected_trait]), ]
  phy_trim <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% traits_trim$species)])
  traits_trim <- traits_trim[match(traits_trim$species, phy_trim$tip.label), ]
  rownames(traits_trim) <- traits_trim$species
  
  # run binaryPGLMM on selected trait
  model <- binaryPGLMM(formula(paste(selected_trait, "~ habit", sep="")), phy=phy_trim, data=traits_trim)
  
  # get model results
  list(trait = selected_trait,
       sigmasq = model$s2, 
       sigmap = model$P.H0.s2, 
       coeff = model$B[2,1], 
       se=model$B.se[2,1], 
       zscore = model$B.zscore[2,1], 
       pval = model$B.pvalue[2,1])
}

# set up input for mapping run_bPGLMM() to list of selected traits
map_input <- list(
  as.list(c("gemmae", "glands", "hairs", "morphotype")),
  list(traits),
  list(phy)
)

# map run_bPGLMM() to selected traits
bpglmm.results <- pmap_dfr(map_input, run_bPGLMM)

# rename rows for xtable
bpglmm.results <- as.data.frame(bpglmm.results)
rownames(bpglmm.results) <- bpglmm.results$trait
bpglmm.results$trait <- NULL

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
