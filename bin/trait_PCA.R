# trait_PCA.R

# Run standard and phylogentic flavors of principal components analysis 
# on traits of ferns from Moorea

# load packages
library(mooreaferns)
library(phytools) # phylo PCA
library(FactoMineR) # standard PCA
library(RColorBrewer) # plotting
library(ggplot2) # plotting
library(cowplot) # plotting
library(ggrepel) # plotting

# set working directory
setwd(here::here())

###########################################################
### function to run PCA and output results for plotting ###
###########################################################

run_trait_PCA <- function (traits, phy, analysis) {
  
  rownames(traits) <- traits$species
  
  # keep only quantitative traits with decent sampling
  sporo_traits <- colnames(traits)[which(sapply(traits, class) == "numeric")]
  
  traits <- traits[,c("habit", sporo_traits)]
  
  # only keep completely sampled data
  traits <- traits[complete.cases(traits), ]
  
  # trim to only species with trait data
  phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])
  
  # get traits in same order as tips
  traits <- traits[phy$tip.label,]
  
  ##### run PCA using FactoMineR package ####
  if (analysis == "standard_PCA") {
    pca.results <- PCA(traits[2:ncol(traits)], graph=FALSE)
    axis1.percent <- round(pca.results$eig[,2][1], digits=1)
    axis2.percent <- round(pca.results$eig[,2][2], digits=1)
    
    # summary of correlaiton between variables and major axes
    # dimdesc(pca.results, axes=c(1,2))
    # summary(pca.results)
    
    # get first two PCA axes
    traits.pca <- as.data.frame(pca.results$ind$coord[,1:2])
    variables.pca <- as.data.frame(pca.results$var$coord[,1:2])
    
    # get contribution to variance of first two PCs
    variance.pca <- as.data.frame(t(pca.results$eig))[2:3,1:2]
    colnames(variance.pca) <- c("Dim.1", "Dim.2")
    rownames(variance.pca) <- c("total_variance", "cumulative_variance")
  }
  
  ###########################################
  ### run phylogenetic PCA using phytools ###
  ###########################################
  if (analysis == "phylo_PCA") {
    phy.pca.results <- phyl.pca(phy, traits[2:ncol(traits)], method="BM")
    phy.pca.summary <- summary(phy.pca.results)
    
    traits.pca <- as.data.frame(phy.pca.results$S[,1:2])
    variables.pca <- as.data.frame(phy.pca.results$L[,1:2])
    variance.pca <- as.data.frame(phy.pca.summary$importance)[2:3,1:2]
    axis1.percent <- round(phy.pca.summary$importance[2,1], digits=4)*100
    axis2.percent <- round(phy.pca.summary$importance[2,2], digits=4)*100
    
    colnames(traits.pca) <- c("Dim.1", "Dim.2")
    colnames(variables.pca) <- c("Dim.1", "Dim.2")
    colnames(variance.pca) <- c("Dim.1", "Dim.2")
    
    rownames(variance.pca) <- c("total_variance", "cumulative_variance")
    variance.pca$Dim.1 <- variance.pca$Dim.1*100
    variance.pca$Dim.2 <- variance.pca$Dim.2*100
  }
  
  traits.pca$species <- rownames(traits.pca)
  traits$species <- rownames(traits)
  traits.pca <- merge(traits.pca, traits[,c("habit", "species")], by="species")
  
  output <- list(traits.pca, variables.pca, axis1.percent, axis2.percent, variance.pca)
  return(output)
}

#######################
### set input files ###
#######################

# load tree
phy <- mooreaferns::fern_tree

# load traits
traits <- mooreaferns::fern_traits

# transform traits
traits <- transform_traits(traits)

###############
### run PCA ###
###############

std_PCA_output <- run_trait_PCA (traits, phy, "standard_PCA")
phy_PCA_output <- run_trait_PCA (traits, phy, "phylo_PCA")

##################################################
### summarize PCA results table for supplement ###
##################################################

std_PCA_loadings <- (std_PCA_output[[2]])
colnames(std_PCA_loadings) <- paste("Std_PCA", colnames(std_PCA_loadings), sep="_")
std_PCA_variance <- (std_PCA_output[[5]])

phy_PCA_loadings <-  (phy_PCA_output[[2]])
colnames(phy_PCA_loadings) <- paste("Phy_PCA", colnames(phy_PCA_loadings), sep="_")
phy_PCA_variance <- (phy_PCA_output[[5]])

PCA.results <- merge(std_PCA_loadings, phy_PCA_loadings, by="row.names")
rownames(PCA.results) <- PCA.results$Row.names
PCA.results$Row.names <- NULL

PCA_variance <- cbind(std_PCA_variance, phy_PCA_variance)
colnames(PCA_variance) <- colnames(PCA.results)

PCA.results <- rbind(PCA.results, PCA_variance)

###########################
### plotting functions  ###
###########################

# function to abbreviate trait names to two or three characters for plotting
abbreviate_traits <- function (trait_name) {
  if (nchar(trait_name) > 3) {
    trait_name <- paste(toupper(substr(trait_name, 1, 1)), substr(trait_name, 2, 2), sep="")
  } else {
    trait_name <- toupper(trait_name)
  }
}

# function to make PCA scores plot
make_PCA_plot <- function (traits.pca, axis1.percent, axis2.percent) {
  p <- ggplot(traits.pca, aes(x=Dim.1, y=Dim.2, color=habit)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_point(aes(colour=habit),size=3,shape=16) +
    labs(x= paste("PC1 (", axis1.percent, "%)", sep=""), y=paste("PC2 (", axis2.percent, "%)", sep="")) +
    standard_theme +
    scale_color_manual(values = cols) +
    theme(legend.position="none")
  return(p)
}

# function to make PCA loadings plot
make_var_plot <- function (variables.pca, axis1.percent, axis2.percent) {
  p <- ggplot(variables.pca, aes(x=Dim.1, y=Dim.2)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_point(,size=3,shape=16) +
    labs(x= paste("PC1 (", axis1.percent, "%)", sep=""), y=paste("PC2 (", axis2.percent, "%)", sep="")) +
    standard_theme +
    geom_text_repel(aes(label = sapply(rownames(variables.pca), abbreviate_traits)), show.legend=FALSE, segment.colour = NA) +
    theme(legend.position="none")
  return(p)
}

##################
### make plots ###
##################

# set colors (epiphytic green, terrestrial brown)
cols <- brewer.pal(9, "Set1")[c(7,3)]

# make standard PCA scores plot
standard_PCA_scores_plot <- make_PCA_plot (std_PCA_output[[1]], std_PCA_output[[3]], std_PCA_output[[4]])

# make standard PCA loadings plot
standard_PCA_loadings_plot <- make_var_plot(std_PCA_output[[2]], std_PCA_output[[3]], std_PCA_output[[4]])
# recenter plot on 0,0
standard_PCA_loadings_plot <- standard_PCA_loadings_plot  +
  expand_limits(y=c(-0.8, 0.8), x=c(-1,1))

# make phylogenetic PCA scores plot
phy_PCA_scores_plot <- make_PCA_plot (phy_PCA_output[[1]], phy_PCA_output[[3]], phy_PCA_output[[4]])
# recenter plot on 0,0
phy_PCA_scores_plot <- phy_PCA_scores_plot +
  expand_limits(y=c(-1, 1), x=c(-1,1))

# make phylogenetic PCA loadings plot
phy_PCA_loadings_plot <- make_var_plot(phy_PCA_output[[2]], phy_PCA_output[[3]], phy_PCA_output[[4]])
# recenter plot on 0,0
phy_PCA_loadings_plot <- phy_PCA_loadings_plot +
  expand_limits(y=c(-0.6, 0.6), x=c(-1,1))

#####################
### combine plots ###
#####################

PCA_plot <- plot_grid(standard_PCA_loadings_plot, 
                      phy_PCA_loadings_plot, 
                      standard_PCA_scores_plot, 
                      phy_PCA_scores_plot, 
                      labels = c("(a)", "(b)", "(c)", "(d)"), 
                      ncol = 2, align="hv")

save_plot("PCA.pdf", PCA_plot, base_height = 5.5, base_width = 6)

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
