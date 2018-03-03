# run Pagel's (1994) test of correlated evolution on two binary traits using fitPagel in phytools

library(xlsx)
library(phytools)

setwd("/Users/joelnitta/R/moorea/")
source("bin/my_funcs.R")

######################
### data wranglin' ###
######################

### get raw data
# set input files
tree_file <- set_input_files()[[2]]

# load tree
phy <- read.tree(tree_file)

# load untransformed traits
traits <- get.my.traits(log.trans = FALSE, scale.traits = FALSE)

# load distribution data (widespread or not)
# distribution data is in "results.df"
source("bin/alpha_diversity/elevational range table.R")

### merge distribution data with traits
# first make all not widespread
traits$distribution <- "notwidespread"
# get list of widespread species
widespread_species <- rownames(results.df[results.df$widespread == 1,])
# make those in list widespread
traits$distribution[rownames(traits) %in% widespread_species] <- "widespread"
rm(results.df)

### match up traits and tree
# trim tree to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])
# get traits in same order as tips
traits <- traits[phy$tip.label,]

##################################
### format traits for fitPagel ###
##################################

# each trait we want to test needs to be a binary, named character vector
traits$habit[traits$habit == 1] <- "epi"
traits$habit[traits$habit == 0] <- "terr"

# make binary morph category: 0 is noncordate, 1 is cordate
traits$morph <- "noncordate"
traits$morph[which(traits$morphotype == "cordate")] <- "cordate"
traits$morphotype <- NULL

traits$glands[traits$glands == 1] <- "glands"
traits$glands[traits$glands == 0] <- "noglands"

traits$hairs[traits$hairs == 1] <- "hairs"
traits$hairs[traits$hairs == 0] <- "nohairs"

traits$gemmae[traits$gemmae == 1] <- "gemmae"
traits$gemmae[traits$gemmae == 0] <- "nogemmae"

# keep only traits of interest
traits.habit <- traits[,c("habit", "morph", "glands", "hairs", "gemmae")]
traits.dist <- traits[,c("distribution", "morph", "glands", "hairs", "gemmae")]

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
data.list.dist <- make_data_list (traits.dist)

####################
### run fitPagel ###
####################

# wrapper to loop through data.list with fitPagel
run_fitPagel <- function (data.list) {
result.list <- list()
for (i in 1:length(data.list)) {
  result <- fitPagel(tree=data.list[[i]]$tree, x=data.list[[i]]$trait, y=data.list[[i]]$indep_var)
  result.list[[i]] <- c(data.list[[i]]$trait.name, signif(result$independent.logL, digits=3), signif(result$dependent.logL, digits=3), signif(result$lik.ratio, digits=3), signif(result$P, digits=3))
}
return (result.list)
}

# run with growth habit as independent var
results <- run_fitPagel (data.list.habit)
summary <-as.data.frame(t(as.data.frame(results)))
rownames(summary) <- NULL
colnames(summary) <- c("trait", "logL_indep", "logL_dep", "likelihood_ratio", "pval")
summary$indep_var <- "habit"

# run with distribution type as independent var
results.dist <- run_fitPagel (data.list.dist)
summary.dist <-as.data.frame(t(as.data.frame(results.dist)))
rownames(summary.dist) <- NULL
colnames(summary.dist) <- c("trait", "logL_indep", "logL_dep", "likelihood_ratio", "pval")
summary.dist$indep_var <- "distribution"

# combine the two
summary <- rbind(summary, summary.dist)

# write out results
write.csv(summary, file = make_filename("correlated_evo_gameto", ".csv", date=TRUE))