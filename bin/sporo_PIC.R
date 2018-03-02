# library(vegan)
# library(picante)
#library(ade4)
library(caper)
library(xlsx)
library(phytools)

setwd("/Users/joelnitta/R/moorea/")
source("bin/my_funcs.R")

# Use brunch to run PIC against binary state (epiphytic or not)

# set input files
tree_file <- set_input_files()[[2]]

# load tree
phy <- read.tree(tree_file)

# load traits
source("bin/traits/get_traits.R")

# keep only quantitative traits with decent sampling
sporo_traits <-c("habit", "stipe", "length", "width", "rhizome", "dissection", "pinna", "SLA")

# keep only quantitative traits with decent sampling
traits <- traits[,sporo_traits]

# make habit a factor
traits$habit[traits$habit == 0] <- "terrestrial"
traits$habit[traits$habit == 1] <- "epiphytic"
traits$habit <- as.factor(traits$habit)

# trim to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])

# get traits in same order as tips
traits <- traits[phy$tip.label,]

########################################################
### phylogenetic ANOVA using phylANOVA from phytools ###
########################################################

# only runs one trait at a time, can't handle NAs
# make loop to trim each trait and run phylANOVA
out <- NULL
groups <- NULL
trait.select <- NULL
habit <- NULL
my.tree <- NULL
result <- list()
for (i in 2:ncol(traits)) {
  trait.select <- traits[,i]
  habit <- traits[,1]
  names(trait.select) <- rownames(traits)
  names(habit) <- rownames(traits)
  trait.select <- trait.select[!(is.na(trait.select))]
  habit <- habit[names(trait.select)]
  my.tree <- phy
  my.tree <- drop.tip(my.tree, my.tree$tip.label[!(my.tree$tip.label %in% names(trait.select))])
  trait.select <- trait.select[my.tree$tip.label]
  habit <- habit[my.tree$tip.label]
  out <- phylANOVA(my.tree, habit, trait.select)
  result[[i-1]] <- c(out$F, out$Pf, colnames(traits)[i])
  names(result[[i-1]]) <- c("F", "p", "trait")
}

summary.panova <- as.data.frame(t(as.data.frame(result)))
rownames(summary.panova) <- NULL

write.csv(summary.panova, make_filename("sporo_phylANOVA", ".csv", date=TRUE))


#######################
### brunch in caper ###
#######################

# caper comparative data command combines tree and data
# but, need column of species names
traits$species <- rownames(traits)
fern<-comparative.data(phy = phy, data = traits, names.col = species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
traits$species <- NULL

# make list of traits to test (drop habit)
trait.list <- colnames(traits)
trait.list <- trait.list[2:ncol(traits)]

# set up loop to run brunch on all traits of interest
result <- list()
model <- NULL
table <- NULL
contrasts <- NULL
num_contrasts <- NULL
num_pos_con <- NULL
tval <- NULL
pval <- NULL

for (i in 1:length(trait.list)) {
  model <- brunch(formula=as.formula(paste(trait.list[i], " ~ habit", sep="")), data=fern)
  table <- caic.table(model)
  contrasts <- table[,1]
  num_contrasts <- length(contrasts)
  num_pos_con <- length(contrasts[contrasts > 0])
  tval <- round(as.vector(summary(model)$coefficients[,3]), digits = 4)
  pval <- round(as.vector(summary(model)$coefficients[,4]), digits = 4)
  result[[i]] <- c(trait.list[i], num_contrasts, num_pos_con, tval, pval)
  names(result[[i]]) <- c("trait", "num_contrasts", "positive_contrasts", "t-value", "p-value")
}

summary.PIC <- as.data.frame(t(as.data.frame(result)))
rownames(summary.PIC) <- NULL

write.csv(summary.PIC, make_filename("sporo_PIC", ".csv", date=TRUE))

# optional: plot the contrasts
# brunchTab <- caic.table(model)
# plot(SLA ~ habit, brunchTab)