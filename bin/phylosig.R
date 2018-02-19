# phylosig.R
# Analyze phylogenetic signal in traits
# Use Blomberg's K and Pagel's lambda for quantitative traits
# and Fritz and Purvis' D for qualitative (binary) traits

# load packages
library(phytools) #phylosig
library(caper) #phylo.d
library(plyr) # relevel
library(mooreaferns)

# load packages
setwd(here::here())

################################################
### load, clean data for quantitative traits ###
################################################

# load tree
phy <- mooreaferns::fern_tree

# load un-transformed species traits
traits <- mooreaferns::fern_traits

# keep only quantitative traits with decent sampling
rownames(traits) <- traits$species
sporo_traits <-c("stipe", "length", "width", "rhizome", "dissection", "pinna", "sla")
traits <- traits[,sporo_traits]

# trim phylogeny to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])

# get traits in same order as tips
traits <- traits[phy$tip.label,]

###########################################
### run Blomberg's K and Pagel's lambda ###
###########################################

# use phytools phylosig function to calculate K, lambda for each trait and summarize in result table
# phylosig will automatically drop species missing trait data

# Blomberg's K
trait.test <- NULL
trait.name <- NULL
phylosig.out <- NULL
kval <- NULL
pval <- NULL
parameters <- NULL
for (i in 1:ncol(traits)) {
  trait.test <- traits[,i]
  names(trait.test) <- rownames(traits)
  phylosig.out <- phylosig(phy, trait.test, method = "K", test = TRUE)
  kval <- phylosig.out$K
  pval <- phylosig.out$P
  trait.name <- colnames(traits)[i]
  parameters[[i]] <- list(trait=trait.name, K = kval, k.pval = pval)
}

k.summary <- as.data.frame(dplyr::bind_rows(parameters))


# Pagel's lambda
trait.test <- NULL
trait.name <- NULL
phylosig.out <- NULL
lambda <- NULL
pval <- NULL
parameters <- NULL
for (i in 1:ncol(traits)) {
  trait.test <- traits[,i]
  names(trait.test) <- rownames(traits)
  phylosig.out <- phylosig(phy, trait.test, method = "lambda", test = TRUE)
  lambda <- phylosig.out$lambda
  pval <- phylosig.out$P
  trait.name <- colnames(traits)[i]
  parameters[[i]] <- list(trait=trait.name, lambda = lambda, lambda.pval = pval)
}

lambda.summary <- as.data.frame(dplyr::bind_rows(parameters))

########################################################
### load, clean data for qualitative (binary) traits ###
########################################################

# starting over, so remove tree and traits for clarity
rm(phy)
rm(traits)

# load tree
phy <- mooreaferns::fern_tree

# load un-transformed species traits
traits <- mooreaferns::fern_traits

# trim to only species with trait data
rownames(traits) <- traits$species
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])

# get traits in same order as tips
traits <- traits[phy$tip.label,]

# make binary morphotype category: 0 is noncordate, 1 is cordate
traits$morph_binary <- 0
traits$morph_binary[which(traits$morphotype == "cordate")] <- 1

##############################
### run Fitz and Purvis' D ###
##############################

# for F&P D value, need to make into caper comparative data object.
gameto.comp <- comparative.data(phy, traits, "species")

# run phylo.d
# would rather do this as a loop but can't figure out how to specify variable
parameters <- list()

phylo.d.out <- phylo.d(gameto.comp, binvar=habit)
parameters[[1]] <- list(phylo.d.out$binvar, phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)
names(parameters[[1]]) <- c("trait", "D", "prob_random", "prob_brownian")

phylo.d.out <- phylo.d(gameto.comp, binvar=glands)
parameters[[2]] <- list(phylo.d.out$binvar, phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)
names(parameters[[2]]) <- c("trait", "D", "prob_random", "prob_brownian")

phylo.d.out <- phylo.d(gameto.comp, binvar=hairs)
parameters[[3]] <- list(phylo.d.out$binvar, phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)
names(parameters[[3]]) <- c("trait", "D", "prob_random", "prob_brownian")

phylo.d.out <- phylo.d(gameto.comp, binvar=gemmae)
parameters[[4]] <- list(phylo.d.out$binvar, phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)
names(parameters[[4]]) <- c("trait", "D", "prob_random", "prob_brownian")

phylo.d.out <- phylo.d(gameto.comp, binvar=morph_binary)
parameters[[5]] <- list(phylo.d.out$binvar, phylo.d.out$DEstimate, phylo.d.out$Pval1, phylo.d.out$Pval0)
names(parameters[[5]]) <- c("trait", "D", "prob_random", "prob_brownian")

d.summary <- as.data.frame(dplyr::bind_rows(parameters))

###########################
### combine the results ###
###########################

phylosig.results <- merge(lambda.summary, k.summary, by="trait")
phylosig.results <- merge(phylosig.results, d.summary, by="trait", all=TRUE)

# rename rows as traits and reorder
phylosig.results$trait <- factor(phylosig.results$trait, levels = c("habit", "stipe", "length", "width", "dissection", "pinna", "sla", "rhizome", "gemmae", "glands", "hairs", "morph_binary") )
rownames(phylosig.results) <- phylosig.results$trait
phylosig.results <- phylosig.results[levels(phylosig.results$trait), ]
phylosig.results$trait <- NULL

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
