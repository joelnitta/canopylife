# microclimate.R
# run ANCOVA comparing epiphytic and terrestrial RH, temp on Moorea
# plot grand mean values with modeled fit

# load packages
library(RColorBrewer) # brewer.pal
library(mooreaferns) # datasets and custom functions
library(tidyverse)
library(cowplot) # plot_grid

# set working directory
setwd(here::here())

# calculate mean microclimate values --------------------------------------

source("bin/microclimate_means.R")

# load and wrangle data ---------------------------------------------------

# load data for sites on Moorea
moorea_sites <- mooreaferns::sites[grep("Aorai", mooreaferns::sites$site, invert=TRUE), ]

# add latitude, longitude, and elevation to microclimate means
daily_values <- left_join(daily_values, moorea_sites)
grand_means <- left_join(grand_means, moorea_sites)

# convert growth habit to factor
daily_values$habit <- factor(daily_values$habit, levels=c("terrestrial", "epiphytic"))
grand_means$habit <- factor(grand_means$habit, levels=c("terrestrial", "epiphytic"))

# run ANCOVA --------------------------------------------------------------

# make function to run ancova on a response variable of choice
run_ancova <- function (resp_var) {
  
  # test for different intercepts between epiphytic and terrestrial dataloggers (but constant slope)
  model1 <- summary(lm(formula(paste(resp_var, " ~ habit + el", sep="")), data=daily_values))
  
  # test for different slopes (and different intercepts) between epiphytic and terrestrial dataloggers 
  model2 <- summary(lm(formula(paste(resp_var, " ~ habit * el", sep="")), data=daily_values))
  
  # compile results
  list (
    resp_var = resp_var,
    
    # pval for different intercepts between epiphytic and terrestrial dataloggers (with same slope)
    pval.intercept = model1$coefficients[2,4],
    f.intercept = model1$fstatistic[1],
    df.intercept = paste(model1$fstatistic[2], model1$fstatistic[3], sep=", "),
    
    # pval for different intercepts between epiphytic and terrestrial (with different slopes)
    pval.slope = model2$coefficients[4,4],
    f.slope = model2$fstatistic[1],
    df.slope = paste(model2$fstatistic[2], model2$fstatistic[3], sep=", ")
    )
}

# set response variables
resp_vars <- c("max_temp", "mean_temp", "min_temp", "sd_temp", "max_RH", "mean_RH", "min_RH", "sd_RH")

# map function to variables
ancova.results <- map_df(resp_vars, run_ancova)

# add rownames
ancova.results <- as.data.frame(ancova.results)
rownames(ancova.results) <- ancova.results$resp_var
ancova.results$resp_var <- NULL

# Get slopes and intercepts for plotting ----------------------------------

# Note that this calculates interaction effects for all variables
# (even though min temp isn't significant)

get_plot_values <- function (resp_var) {
  
  model1 <- lm(formula(paste(resp_var, " ~ habit/el + 0", sep="")), data=daily_values)
  model2 <- lm(formula(paste(resp_var, " ~ habit*el", sep="")), data=daily_values)
  
  list(
  resp_var = resp_var,
  slope.epi = coef(model1)[4],
  inter.epi = coef(model1)[2],
  slope.ter = coef(model1)[3],
  inter.ter = coef(model1)[1],
  r.squared = summary(model2)$r.squared,
  pval = summary(model2)$coefficients[4,4] )
  
}

# map function to variables
slopes.results <- map_df(resp_vars, get_plot_values)

# add rownames
slopes.results <- as.data.frame(slopes.results)
rownames(slopes.results) <- slopes.results$resp_var
slopes.results$resp_var <- NULL

# Plotting ----------------------------------------------------------------

# function to make scatterplot with lines and R2 values from linear model
make_scatter_plot <- function (input.plot.data, model_data, indep_var, dep_var, ylabel) {
  
  # extract formatted r squared value for this plot
  r.squared <- paste("italic(R)^2 ==", round(model_data[dep_var,"r.squared"], 2), sep=" " )    
  
  p <- ggplot (data = input.plot.data, aes_string(x = indep_var, y = dep_var, color = "habit")) +
    geom_abline(intercept = model_data[dep_var,"inter.ter"], slope = model_data[dep_var,"slope.ter"], color=cols[1]) +
    geom_abline(intercept = model_data[dep_var,"inter.epi"], slope = model_data[dep_var,"slope.epi"], color=cols[2]) +
    annotate("text",
             x=Inf,
             y=Inf,
             label=r.squared,
             parse=TRUE,
             hjust=1,
             vjust=1) +
    geom_point(size = 2) +
    standard_theme +
    theme(legend.position = "none") +
    xlab("Elevation (m)") +
    ylab(ylabel) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks=c(0,250,500,750,1000))
  return(p)
}

# set colors for manual line plotting: green for epiphytes, brown for terrestrial
qualcols <- brewer.pal(9, "Set1")
cols <- qualcols[c(7,3)]

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

# clean up workspace (reserve "result" in object name only for final results objects to keep)
rm(list=ls()[grep("result", ls(), invert=TRUE)])
