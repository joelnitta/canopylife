# binaryPGLMM.R
# Construct binary phylogenetic generalized linear mixed models (PGLMMs) 
# for binary (gametophyte) traits related to epiphytic growth

# load packages
library(plyr)
library(ape)
library(mooreaferns)

# set working directory
setwd(here::here())

# clear workspace
rm(list=ls())

########################################################
### load, clean data for qualitative (binary) traits ###
########################################################

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

####################
### binaryPGLMM ####
####################

# define traits to test
traits.test <- c("glands", "hairs", "gemmae", "morph_binary")

# binaryPGLMM can't take NA values
# make list of trait dataframes and phylogenies which only include species with no NA values for each of traits to test
traits.df <- list()
phy.list <- list()
for (i in 1:length(traits.test)) {
  traits.df[[i]] <- traits[!(is.na(traits[,traits.test[i]])), ]
  phy.list[[i]] <- drop.tip(phy,  phy$tip.label[!(phy$tip.label %in% rownames(traits.df[[i]]))])
  traits.df[[i]] <- traits.df[[i]][phy.list[[i]]$tip.label, ]
}

# run binary PGLMM as loop
model <- list()
results <- list()
for (i in 1:length(traits.test)) {
  model[[i]] <- binaryPGLMM(formula(paste(traits.test[i], "~ habit", sep="")), phy=phy.list[[i]], data=traits.df[[i]] )
  results[[i]] <- list(sigmasq = model[[i]]$s2, 
                       sigmap = model[[i]]$P.H0.s2, 
                       coeff = model[[i]]$B[2,1], 
                       se=model[[i]]$B.se[2,1], 
                       zscore = model[[i]]$B.zscore[2,1], 
                       pvalue = model[[i]]$B.pvalue[2,1])
  }

# compile results into single table
bpglmm.summary <- do.call(rbind.data.frame, results)

### format for manuscipt
# round to 3 digits
bpglmm.summary <- as.data.frame(sapply(bpglmm.summary, round, 3))

# rename rows as traits and reorder alphabetically
bpglmm.summary$trait <- traits.test
bpglmm.summary$trait <- factor(bpglmm.summary$trait, levels = sort(bpglmm.summary$trait))
bpglmm.summary$trait <- mapvalues(bpglmm.summary$trait, 
                                    c("glands", "hairs", "gemmae", "morph_binary"),
                                    c("Glands", "Hairs", "Gemmae", "Morphotype"))
rownames(bpglmm.summary) <- bpglmm.summary$trait
bpglmm.summary <- bpglmm.summary[levels(bpglmm.summary$trait), ]
bpglmm.summary$trait <- NULL
