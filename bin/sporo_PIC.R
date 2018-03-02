# sporoPIC

# Calculate phylogenetically independent contrasts for quantitative (sporophyte) traits 
# between terrestrial and epiphytic species using brunch function in caper

# load packages
library(tidyverse)
library(caper)
library(mooreaferns)

# set working directory
setwd(here::here())

####################
### prepare data ###
####################

# load tree
phy <- mooreaferns::fern_tree

# load traits
traits <- mooreaferns::fern_traits

# reset levels for growth habit so epiphytic trait values are subtracted from terrestrial
traits$habit <- factor(traits$habit, levels = c("terrestrial", "epiphytic"))

# keep only species, habit, and quantitative traits to test
traits.test <- c("stipe", "length", "width", "rhizome", "dissection", "pinna", "sla")
traits <- traits[,c("species", "habit", traits.test)]
rownames(traits) <- traits$species

# trim to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])

# get traits in same order as tips
traits <- traits[phy$tip.label,]

##############################
### run PICs using brunch  ###
##############################

### brunch requires caper comparative data objects (combining tree/trait datasets) as input 
# I'm not sure how NA values affect brunch calculation of PICs, so make list of caper comparative 
# data objects, each one with NAs in traits removed
fern_data_list <- list()
for (i in 1:length(traits.test)) {
  fern_data_list[[i]] <- comparative.data(phy, traits[,c("species", "habit", traits.test[i])], "species", na.omit=TRUE)
}

# set up loop to run brunch on each trait of interest, one at a time
result <- list()
model <- NULL
table <- NULL
contrasts <- NULL
num_contrasts <- NULL
num_pos_con <- NULL
tval <- NULL
pval <- NULL

# run loop
for (i in 1:length(fern_data_list)) {
  model <- brunch(formula=as.formula(paste(traits.test[i], " ~ habit", sep="")), data=fern_data_list[[i]])
  table <- caic.table(model)
  contrasts <- table[,1]
  num_contrasts <- length(contrasts)
  num_pos_con <- length(contrasts[contrasts > 0])
  tval <- as.vector(summary(model)$coefficients[,3])
  pval <- as.vector(summary(model)$coefficients[,4])
  result[[i]] <- c(num_contrasts, num_pos_con, tval, pval)
  names(result[[i]]) <- c("num_contrasts", "num_pos_con", "tval", "pval")
}
names(result) <- traits.test

# combine results
PIC.results <- bind_rows(result)

# transpose
PIC.results %>%
  rownames_to_column %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) ->
  PIC.results

# store as dataframe (need rownames for xtable)
PIC.results <- as.data.frame(PIC.results)

# format results dataframe
rownames(PIC.results) <- PIC.results$var
PIC.results$var <- NULL
colnames(PIC.results) <- c("num_contrasts", "num_pos_con", "tval", "pval")
PIC.results <- PIC.results[c("stipe", "length", "width", "dissection", "pinna", "sla", "rhizome"), ]
PIC.results$num_contrasts <- as.integer(PIC.results$num_contrasts)
PIC.results$num_pos_con <- as.integer(PIC.results$num_pos_con)
