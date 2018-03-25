# cont_phylosig.R

# Analyze phylogenetic signal in continuous traits
# using Blomberg's K and Pagel's lambda

# load packages
library(dplyr) # bind_rows
library(phytools) # phylosig
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

k.summary <- as.data.frame(bind_rows(parameters))


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

lambda.summary <- as.data.frame(bind_rows(parameters))

###########################
### combine the results ###
###########################

cont_phylosig.results <- merge(lambda.summary, k.summary, by="trait")

# name rows by trait and reorder
rownames(cont_phylosig.results) <- cont_phylosig.results$trait
cont_phylosig.results$trait <- NULL
cont_phylosig.results <- cont_phylosig.results[c("stipe", "length", "width", "dissection", "pinna", "sla", "rhizome"), ]

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
