# binary_phylosig.R

# Analyze phylogenetic signal in binary (qualitative) traits
# using Fritz and Purvis' D

# load packages
library(caper) # phylo.d()
library(mooreaferns)

# set working directory
setwd(here::here())

###########################
### load and clean data ###
###########################

# load tree
phy <- mooreaferns::fern_tree

# load un-transformed species traits
traits <- mooreaferns::fern_traits

# make cordate morphotype category: 0 is noncordate, 1 is cordate
traits$cordate_morph <- 0
traits$cordate_morph[which(traits$morphotype == "cordate")] <- 1

# make growth habit into a binary category: 0 is not epiphtyic (ie, terrestrial), 1 is epiphytic
# growth habit is already a 2-level factor, so just convert to numeric and adjust
traits$habit <- as.numeric(traits$habit)
traits$habit <- traits$habit-1

### make caper tree/trait datasets silently dropping missing data
# some of the traits include NAs.
# for phylo.d, need to have matching tree and trait data with all NAs removed
# don't do comparative.data on all the traits together, or species missing data for ANY trait will get dropped
habit.comp  <- comparative.data(phy, traits[,c("species", "habit")], "species", na.omit=TRUE)
gemmae.comp <- comparative.data(phy, traits[,c("species", "gemmae")], "species", na.omit=TRUE)
glands.comp <- comparative.data(phy, traits[,c("species", "glands")], "species", na.omit=TRUE)
hairs.comp  <- comparative.data(phy, traits[,c("species", "hairs")], "species", na.omit=TRUE)
morph.comp  <- comparative.data(phy, traits[,c("species", "cordate_morph")], "species", na.omit=TRUE)

##############################
### run Fitz and Purvis' D ###
##############################

# I would rather do this as a loop but can't figure out how to specify binvar, so hard-code it

# make list to store observed D and associated P values for each trait
parameters <- list()

# run phylo.d on a single trait, store output as first item in list
phylo.d.out <- phylo.d(habit.comp, binvar=habit)
parameters[[1]] <- list(phylo.d.out$binvar, phylo.d.out$StatesTable[1], phylo.d.out$StatesTable[2], phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)

# repeat for other traits
phylo.d.out <- phylo.d(gemmae.comp, binvar=gemmae)
parameters[[2]] <- list(phylo.d.out$binvar, phylo.d.out$StatesTable[1], phylo.d.out$StatesTable[2], phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)

phylo.d.out <- phylo.d(glands.comp, binvar=glands)
parameters[[3]] <- list(phylo.d.out$binvar, phylo.d.out$StatesTable[1], phylo.d.out$StatesTable[2], phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)

phylo.d.out <- phylo.d(hairs.comp, binvar=hairs)
parameters[[4]] <- list(phylo.d.out$binvar, phylo.d.out$StatesTable[1], phylo.d.out$StatesTable[2], phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)

phylo.d.out <- phylo.d(morph.comp, binvar=cordate_morph)
parameters[[5]] <- list(phylo.d.out$binvar, phylo.d.out$StatesTable[1], phylo.d.out$StatesTable[2], phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)

# combine results into single dataframe
parameters <- lapply(parameters, setNames, c("trait", "num_present", "num_absent", "D", "prob_random", "prob_brownian"))
binary_phylosig.results <- as.data.frame(dplyr::bind_rows(parameters))

# name rows by trait and reorder
rownames(binary_phylosig.results) <- binary_phylosig.results$trait
binary_phylosig.results$trait <- NULL
binary_phylosig.results <- binary_phylosig.results[c("habit", "gemmae", "glands", "hairs"), ]

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
