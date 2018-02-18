# Analyze phylogenetic signal in traits
# Use Blomberg's K and Pagel's lambda for quantitative traits
# and Fritz and Purvis' D for qualitative (binary) traits

library(phytools) #phylosig
library(caper) #phylo.d
library(plyr) # relevel
library(mooreaferns)

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

# k
trait.test <- NULL
result <- NULL
k.row.result <- NULL
pval <- NULL
asterisk <- NULL
for (i in 1:ncol(traits)) {
  trait.test <- traits[,i]
  names(trait.test) <- rownames(traits)
  result[[i]] <- phylosig(phy, trait.test, method = "K", test = TRUE)
  pval[i] <- result[[i]]$P
  
  if (is.na(pval[i])) {
    asterisk[i] <- ""
  } else if (pval[i] > 0 & pval[i] < 0.0005) {
    asterisk[i] <- "***"
  } else if (pval[i] >= 0.0005 & pval[i] < 0.005) {
    asterisk[i] <- "**"
  } else if (pval[i] >= 0.005 & pval[i] < 0.05) { 
    asterisk[i] <- "*"
  } else if (pval[i] >= 0.05) { 
    asterisk[i] <- " (NS)"
  }
  
  k.row.result[[i]] <- c(colnames(traits)[i], paste (round(result[[i]]$K, digits = 2), asterisk[i], sep=""))
}

k.summary <-as.data.frame(t(as.data.frame(k.row.result)), row.names = NULL)
rownames(k.summary) <- NULL
colnames(k.summary) <- c("trait", "K")

# lambda
trait.test <- NULL
result <- NULL
l.row.result <- NULL
pval <- NULL
asterisk <- NULL
for (i in 1:ncol(traits)) {
  trait.test <- traits[,i]
  names(trait.test) <- rownames(traits)
  result[[i]] <- phylosig(phy, trait.test, method = "lambda", test = TRUE)
  pval[i] <- result[[i]]$P
  if (is.na(pval[i])) {
    asterisk[i] <- ""
  } else if (pval[i] > 0 & pval[i] < 0.0005) {
    asterisk[i] <- "***"
  } else if (pval[i] >= 0.0005 & pval[i] < 0.005) {
    asterisk[i] <- "**"
  } else if (pval[i] >= 0.005 & pval[i] < 0.05) { 
    asterisk[i] <- "*"
  } else if (pval[i] >= 0.05) { 
    asterisk[i] <- " (NS)"
  }
  l.row.result[[i]] <- c(colnames(traits)[i], paste (round(result[[i]]$lambda, digits = 2), asterisk[i], sep=""))
}
l.summary <-as.data.frame(t(as.data.frame(l.row.result)))
rownames(l.summary) <- NULL
colnames(l.summary) <- c("trait", "lambda")

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
result.list <- list()

result <- phylo.d(gameto.comp, binvar=habit)
result.list[[1]] <- c(result$binvar, result$DEstimate, result$Pval1, result$Pval0)
names(result.list[[1]]) <- c("trait", "Estimated_D", "prob_random", "prob_brownian")

result <- phylo.d(gameto.comp, binvar=glands)
result.list[[2]] <- c(result$binvar, result$DEstimate, result$Pval1, result$Pval0)
names(result.list[[2]]) <- c("trait", "Estimated_D", "prob_random", "prob_brownian")

result <- phylo.d(gameto.comp, binvar=hairs)
result.list[[3]] <- c(result$binvar, result$DEstimate, result$Pval1, result$Pval0)
names(result.list[[3]]) <- c("trait", "Estimated_D", "prob_random", "prob_brownian")

result <- phylo.d(gameto.comp, binvar=gemmae)
result.list[[4]] <- c(result$binvar, result$DEstimate, result$Pval1, result$Pval0)
names(result.list[[4]]) <- c("trait", "Estimated_D", "prob_random", "prob_brownian")

result <- phylo.d(gameto.comp, binvar=morph_binary)
result.list[[5]] <- c(result$binvar, result$DEstimate, result$Pval1, result$Pval0)
names(result.list[[5]]) <- c("trait", "Estimated_D", "prob_random", "prob_brownian")

d.summary <- as.data.frame(t(as.data.frame(result.list)))
rownames(d.summary) <- NULL

d.summary$D <- paste(round(as.numeric(as.character(d.summary$Estimated_D)), 2), " (Prnd = ", d.summary$prob_random, ", Pbm = ", d.summary$prob_brownian, ")", sep="")
d.summary <- d.summary[,c("trait", "D")]

###########################
### combine the results ###
###########################

phylosig.summary <- merge(l.summary, k.summary, by="trait")
phylosig.summary <- merge(phylosig.summary, d.summary, by="trait", all=TRUE)

# adjust classes of columns, rename and reorder traits for manuscript
phylosig.summary[c("lambda", "K", "D")] <- sapply(phylosig.summary[c("lambda", "K", "D")], as.character)

phylosig.summary$trait <- factor(phylosig.summary$trait, levels = c("habit", "stipe", "length", "width", "dissection", "pinna", "sla", "rhizome", "gemmae", "glands", "hairs", "morph_binary") )

phylosig.summary$trait <- mapvalues(phylosig.summary$trait, 
                                  c("dissection", "length", "pinna", "rhizome", "sla", "stipe", "width", "gemmae", "glands", "habit", "hairs", "morph_binary"),
                                  c("Frond Dissection", "Frond Length", "Pinna Number", "Rhizome Diam.", "SLA", "Stipe Length", "Frond Width", "Gemmae", "Glands", "Growth Habit", "Hairs", "Morphotype"))

rownames(phylosig.summary) <- phylosig.summary$trait
phylosig.summary <- phylosig.summary[levels(phylosig.summary$trait), ]

phylosig.summary$trait <- NULL

rm(l.summary, k.summary, d.summary, traits)

# uncomment to write out results as csv
# write.csv(phylosig.summary, "table2_phylogenetic_signal.csv")
