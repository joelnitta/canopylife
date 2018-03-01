# run general linear mixed model including phylogeny to test quant traits vs. growth habit

# to install new package on cluster:
# install.packages("MCMCglmm", lib="~/apps/R", repos="http://cran.us.r-project.org")

library(ape)
library(MCMCglmm)

# working directory on laptop
# setwd("/Users/joelnitta/R/moorea/")

# workign directory on cluster
setwd("/n/home04/jnitta/Analysis/moorea")

# define function to extract p value table from summary as dataframe, append model name
summary_to_csv <- function (model) {
  results.df <- as.data.frame(summary(model)$solutions)
  results.df$model <- paste(model$Fixed$formula[2], model$Fixed$formula[1], model$Fixed$formula[3], sep = " ")
  results.df$effect <- rownames(results.df)
  rownames(results.df) <- NULL
  colnames(results.df) <- c("Parameter estimate", "Lower 95% CI", "Upper 95% CI", "Effective sample size", "P-value", "Model", "Effect")
  return(results.df)
}

# set tree file
# tree_file <- set_input_files()[[2]]
tree_file <- "data/treepl_Moorea_Tahiti_2016-01-03.tre"

# load tree
phy <- read.tree(tree_file)

# run script to calculate trait means, combine with qualitative traits
# traits <- get.my.traits()
# to run on cluster, use pre-processed trait data
trait_file <- "data/traits_glmm.csv"
traits <- read.csv(trait_file, row.names=1)

# trim to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])

# get traits in same order as tips
traits <- traits[phy$tip.label,]

###################
### set up data for MCMCglmm 
##################

# define traits to test
sporo_traits <-c("stipe", "length", "width", "rhizome", "dissection", "pinna", "SLA")
traits.test <- c(sporo_traits, "habit")

# MCMCglmm can't take NA values
# make list of trait dataframes and phylogenies which only include species with no NA values for the traits I'm going to test
traits.list <- list()
phy.list <- list()
inv.phy <- list()
for (i in 1:length(traits.test)) {
  traits.list[[i]] <- traits[!(is.na(traits[,traits.test[i]])), ]
  phy.list[[i]] <- drop.tip(phy,  phy$tip.label[!(phy$tip.label %in% rownames(traits.list[[i]]))])
  inv.phy[[i]] <-inverseA(phy.list[[i]],nodes="TIPS",scale=TRUE)
  traits.list[[i]] <- traits.list[[i]][phy.list[[i]]$tip.label, ]
  traits.list[[i]]$species <- rownames(traits.list[[i]])
}

# run MCMCglmm

# set number of generations
# 10,000 thousand takes about 13 sec
# 5,000,000 takes about 2 hr? (recommended number in tutorial)
# http://www.mpcm-evolution.org/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
# set to 600,000 (about 1 hr) to test on cluster
# num <- 600000
num <- 10000

# define new prior list. Need one G for each random effect, one R for each fixed effect
my.prior<-list(G=list(G1=list(V=1,nu=0.02)),
               R=list(R1=list(V=1,nu=0.02)))

model <- list()
ptm <- NULL
results <- list()
time <- list()
summary <- list()
for (i in 1:(length(traits.test)-1) ) {
  ptm <- proc.time()
  model[[i]] <-MCMCglmm(formula(paste(traits.test[i], "~ habit"), sep=""), 
                        random=~species, 
                        family="gaussian", 
                        ginverse=list(species=inv.phy[[i]]$Ainv), 
                        prior=my.prior, 
                        data=traits.list[[i]], 
                        nitt=num,burnin=1000,thin=500)
  # print processing time
  time[[i]] <- (proc.time() - ptm)[3]
  # print summary
  summary[[i]] <- summary(model[[i]])
  # print traces
  pdf(file=paste("bin/",Sys.Date(),"/","mcmcglmm_trace_", traits.test[i], ".pdf",sep=""))
    plot(model[[i]]$Sol, auto.layout=F)
  dev.off()
  # save output as df
  results[[i]] <- summary_to_csv(model[[i]])
}

# compile list of dfs into single df
final.df <- do.call('rbind', results)

# check processing time
time

# print summary of models
summary

# save final results df to CSV
write.csv (final.df, file=paste("bin/",Sys.Date(),"/","mcmcglmm_habit_results.csv",sep=""))