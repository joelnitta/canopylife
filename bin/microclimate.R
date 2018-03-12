# microclimate.R
# run ANCOVA comparing epiphytic and terrestrial RH, temp on Moorea
# plot grand mean values with modeled fit

# load packages
library(ggplot2) # plotting functions
library(cowplot) # plot_grid
library(RColorBrewer) # brewer.pal
library(mooreaferns) # datasets and custom functions

# set working directory
setwd(here::here())

##########################################
### calculate mean microclimate values ###
##########################################

source("bin/microclimate_means.R")

##################
### run ANCOVA ###
##################

# add latitude, longitude, and elevation to microclimate means
daily_values <- merge(daily_values, moorea_sites, by="site", all.x=TRUE)
grand_means <- merge(grand_means, moorea_sites, by="site", all.x=TRUE)

# check to make sure it worked
#anyNA(daily_values)
#anyNA(grand_means)

# make sure growth habit is factor
daily_values$habit <- as.factor(daily_values$habit)
grand_means$habit <- as.factor(grand_means$habit)

# set response variables
resp_vars <- c("max_temp", "mean_temp", "min_temp", "sd_temp", "max_RH", "mean_RH", "min_RH", "sd_RH")

# set up dummy variables for loop
pval.intercept <- NULL
f.intercept <- NULL
df.intercept <- NULL
pval.slope <- NULL
f.slope  <- NULL
df.slope <- NULL
model1 <- NULL
model2 <- NULL
parameters <- NULL

# run the ANCOVA for each response variable
for (i in 1:length(resp_vars)) {
  # test for different intercepts (but constant slope)
  model1 <- summary(lm(formula(paste(resp_vars[i], " ~ habit + el", sep="")), data=daily_values))
  # pval for different intercepts between epi and ter (with same slope)
  pval.intercept <- model1$coefficients[2,4]
  f.intercept <- model1$fstatistic[1]
  df.intercept <- paste(model1$fstatistic[2], model1$fstatistic[3], sep=", ")
  
  # test for different slopes (and different intercepts)
  model2 <- summary(lm(formula(paste(resp_vars[i], " ~ habit * el", sep="")), data=daily_values))
  pval.slope <- model2$coefficients[4,4]
  f.slope <- model2$fstatistic[1]
  df.slope <- paste(model2$fstatistic[2], model2$fstatistic[3], sep=", ")
  
  parameters[[i]] <- list(df.intercept = df.intercept, f.intercept = f.intercept, pval.intercept = pval.intercept, 
                       df.slope = df.slope, f.slope = f.slope, pval.slope = pval.slope)
}

# combine results into single dataframe
ancova.results <- as.data.frame(dplyr::bind_rows(parameters))
rownames(ancova.results) <- resp_vars

# write.csv(ancova.results, file="climate_ancova.csv")

###############################################
#### Get slopes and intercepts for plotting ###
###############################################

# Note that this calculates interaction effects for all variables
# (even though min temp isn't significant)

model1<- NULL
model2 <- NULL
parameters <- NULL
slope.epi <- list()
inter.epi <- list()
slope.ter <- list()
inter.ter<- list()
r.squared <- list()
pval <- list ()
for (i in 1:length(resp_vars)) {
  model1 <- lm(formula(paste(resp_vars[i], " ~ habit/el + 0", sep="")), data=daily_values)
  model2 <- lm(formula(paste(resp_vars[i], " ~ habit*el", sep="")), data=daily_values)
  parameters <- coef(model1)
  slope.epi[[i]] <- parameters[3]
  inter.epi[[i]] <- parameters[1]
  slope.ter[[i]] <- parameters[4]
  inter.ter[[i]] <- parameters[2]
  r.squared[[i]] <- summary(model2)$r.squared
  pval[[i]] <- summary(model2)$coefficients[4,4]
}

slopes.results <- data.frame(r.squared = unlist(r.squared), slope.epi = unlist(slope.epi), inter.epi=unlist(inter.epi), slope.ter=unlist(slope.ter), inter.ter=unlist(inter.ter), pval = unlist(pval))
rownames(slopes.results) <- resp_vars

################
### Plotting ### 
################

# function to make scatterplot with lines and R2 values from linear model
make_scatter_plot <- function (input.plot.data, model_data, indep_var, dep_var, ylabel) {
  p <- ggplot (data = input.plot.data, aes_string(x = indep_var, y = dep_var, color = "habit")) +
    geom_abline(intercept = model_data[dep_var,"inter.ter"], slope = model_data[dep_var,"slope.ter"], color=cols[1]) +
    geom_abline(intercept = model_data[dep_var,"inter.epi"], slope = model_data[dep_var,"slope.epi"], color=cols[2]) +
    annotate("text", x=Inf, y=Inf, hjust=1, vjust=1, label = paste("R2 = ", round(model_data[dep_var,"r.squared"], 2), sep="")) +
    geom_point(size = 2) +
    standard_theme +
    theme(legend.position = "none") +
    xlab("Elevation (m)") +
    ylab(ylabel) +
    scale_colour_manual(values = cols)
  return(p)
}

# set colors for manual line plotting: green for epiphytes, brown for terrestrial
qualcols <- brewer.pal(9, "Set1")
cols <- qualcols[c(7,3)]

# reset growth habit levels
grand_means$habit <- factor(grand_means$habit, levels =c("ter", "epi"))
levels(grand_means$habit) <- c("terrestrial", "epiphytic")

# make each figure part
partA <- make_scatter_plot (grand_means, slopes.results, "el", "mean_temp", "Mean Temp. (C)")
partA <- partA + blank_x.theme

partB <- make_scatter_plot (grand_means, slopes.results, "el", "sd_temp", "SD Temp. (C)")
partB <- partB + blank_x.theme

partC <- make_scatter_plot (grand_means, slopes.results, "el", "min_RH", "Minimum RH (%)")

partD <- make_scatter_plot (grand_means, slopes.results, "el", "sd_RH", "SD RH (%)")

# assemble the figure parts together
climate_plot <- plot_grid(partA, partB, partC, partD, nrow=2, ncol=2, labels = NULL, align = "hv")

# write out pdf
save_plot("microclimate.pdf", climate_plot, base_height = 4, base_width = 6)

# don't need pvalue and rsquared in slopes.results anymore (already in ancova.results or plot)
slopes.results$pval <- NULL
slopes.results$r.squared <- NULL

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
