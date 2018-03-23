# correlated_evo

# run Pagel's (1994) test of correlated evolution on two binary traits using fitPagel in phytools

# load packages
library(phytools)
library(mooreaferns)

#####################
### prepare data  ###
#####################

# load tree
phy <- mooreaferns::fern_tree

# load un-transformed species traits
traits <- mooreaferns::fern_traits
rownames(traits) <- traits$species

### match up traits and tree
# trim tree to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])
# get traits in same order as tips
traits <- traits[phy$tip.label,]

##################################
### format traits for fitPagel ###
##################################

# make binary morph category: 0 is noncordate, 1 is cordate
traits$morph_binary <- "noncordate"
traits$morph_binary[which(traits$morphotype == "cordate")] <- "cordate"
traits$morphotype <- NULL

traits$glands[traits$glands == 1] <- "glands"
traits$glands[traits$glands == 0] <- "noglands"

traits$hairs[traits$hairs == 1] <- "hairs"
traits$hairs[traits$hairs == 0] <- "nohairs"

traits$gemmae[traits$gemmae == 1] <- "gemmae"
traits$gemmae[traits$gemmae == 0] <- "nogemmae"

# keep only traits of interest
traits.habit <- traits[,c("habit", "morph_binary", "glands", "hairs", "gemmae")]

# fitPagel crashes if there are NAs
# so need to make list of trees and traits with NAs trimmed out
# all traits must be named vectors
# each item in list contains matching indepedent var (habit or distribution), dependent trait, and tree
# first trait in trait df must be independent var
make_data_list <- function (traits) {
  data.list <- vector(mode="list", length=ncol(traits)-1)
  for (i in 1:(ncol(traits)-1)) {
    traits.trim <- traits[!(is.na(traits[,i+1])), ]
    data.list[[i]]$indep_var <- traits.trim[,1]
    names(data.list[[i]]$indep_var) <- rownames(traits.trim)
    data.list[[i]]$trait <- traits.trim[,i+1]
    names(data.list[[i]]$trait) <- rownames(traits.trim)
    data.list[[i]]$tree <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits.trim))])
    data.list[[i]]$trait.name <- colnames(traits.trim)[i+1]
  }
  return (data.list)
}

# make two lists, one for habit and one for distribution type
data.list.habit <- make_data_list (traits.habit)

####################
### run fitPagel ###
####################

# wrapper to loop through data.list with fitPagel
run_fitPagel <- function (data.list) {
  model.list <- list()
  for (i in 1:length(data.list)) {
    model <- fitPagel(tree=data.list[[i]]$tree, x=data.list[[i]]$trait, y=data.list[[i]]$indep_var)
    model.list[[i]] <- list(trait = data.list[[i]]$trait.name, 
                            logL_indep = model$independent.logL, 
                            logL_dep = model$dependent.logL, 
                            likelihood_ratio = model$lik.ratio, 
                            pval = model$P)
  }
  return (model.list)
}

# run with growth habit as independent var
results <- run_fitPagel (data.list.habit)

# combine results
pagel.results <- dplyr::bind_rows(results)

# store as dataframe (need rownames for xtable)
pagel.results <- as.data.frame(pagel.results)

# add rownames
rownames(pagel.results) <- pagel.results$trait
pagel.results$trait <- NULL
