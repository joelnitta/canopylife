# microclimate.R
# run ANCOVA comparing epiphytic and terrestrial RH, temp on Moorea
# plot grand mean values with modeled fit

# load packages
library(ggplot2) # plotting functions
library(cowplot) # plot_grid
library(RColorBrewer) # brewer.pal
library(mooreaferns) # datasets and custom functions
library(plyr) #ddply

# set working directory
setwd(here::here())

# clear workspace
rm(list=ls())

################
### Get Data ###
################

# load site data, using only sites on Moorea
moorea_sites <- mooreaferns::sites[grep("Aorai", sites$site, invert=TRUE), ]

# reset row names
rownames(moorea_sites) <- NULL

# load raw microclimate data for Moorea (temperature and RH every 15 minutes)
climate.data <- (mooreaferns::moorea_climate)

# calculate daily maximum, mean, minimum, and SD of temperature and RH by site
daily_values <- ddply(climate.data, .(site, date), summarize,
                     max_temp = max(temp),
                     mean_temp = mean(temp),
                     min_temp = min(temp),
                     sd_temp = sd(temp),
                     max_RH = max(RH),
                     mean_RH = mean(RH),
                     min_RH = min(RH),
                     sd_RH = sd(RH))

# calculate grand means of daily values for temperature and RH by site
grand_means <- ddply(daily_values, "site", summarize,
                     max_temp = mean(max_temp),
                     mean_temp = mean(mean_temp),
                     min_temp = mean(min_temp),
                     sd_temp = mean(sd_temp),
                     max_RH = mean(max_RH),
                     mean_RH = mean(mean_RH),
                     min_RH = mean(min_RH),
                     sd_RH = mean(sd_RH))

### reformat data to be in site x variable format
# add a column for growth habit: daily means
daily_values$habit <- ""
daily_values$habit[grep("epi", daily_values$site)] <- "epi"
daily_values$habit[grep("ter", daily_values$site)] <- "ter"
daily_values$site <- gsub("_epi", "", daily_values$site)
daily_values$site <- gsub("_ter", "", daily_values$site)
# grand means
grand_means$habit <- ""
grand_means$habit[grep("epi", grand_means$site)] <- "epi"
grand_means$habit[grep("ter", grand_means$site)] <- "ter"
grand_means$site <- gsub("_epi", "", grand_means$site)
grand_means$site <- gsub("_ter", "", grand_means$site)

# add latitude, longitude, and elevation
daily_values <- merge(daily_values, moorea_sites, by="site", all.x=TRUE)
grand_means <- merge(grand_means, moorea_sites, by="site", all.x=TRUE)

# check to make sure it worked
#anyNA(daily_values)
#anyNA(grand_means)

# make sure growth habit is factor
daily_values$habit <- as.factor(daily_values$habit)
grand_means$habit <- as.factor(grand_means$habit)

rm(climate.data)

##################
### run ANCOVA ###
##################

# set response variables
resp_vars <- c("max_temp", "mean_temp", "min_temp", "sd_temp", "max_RH", "mean_RH", "min_RH", "sd_RH")

# set up dummy variables for loop
pval.intercept <- NULL
f.intercept <- NULL
df.intercept <- NULL
pval.slope <- NULL
f.slope  <- NULL
df.slope <- NULL
out1 <- NULL
out2 <- NULL
results <- NULL

# run the ANCOVA for each response variable
for (i in 1:length(resp_vars)) {
  # test for different intercepts (but constant slope)
  out1 <- summary(lm(formula(paste(resp_vars[i], " ~ habit + el", sep="")), data=daily_values))
  # pval for different intercepts between epi and ter (with same slope)
  pval.intercept <- out1$coefficients[2,4]
  f.intercept <- out1$fstatistic[1]
  df.intercept <- paste(out1$fstatistic[2], out1$fstatistic[3], sep=", ")
  
  # test for different slopes (and different intercepts)
  out2 <- summary(lm(formula(paste(resp_vars[i], " ~ habit * el", sep="")), data=daily_values))
  pval.slope <- out2$coefficients[4,4]
  f.slope <- out2$fstatistic[1]
  df.slope <- paste(out2$fstatistic[2], out2$fstatistic[3], sep=", ")
  
  results[[i]] <- list(df.intercept = df.intercept, f.intercept = f.intercept, pval.intercept = pval.intercept, 
                       df.slope = df.slope, f.slope = f.slope, pval.slope = pval.slope)
}

# combine results into single dataframe
ancova_results <- do.call(rbind.data.frame, results)
rownames(ancova_results) <- resp_vars

# write.csv(results, file=make_filename("climate_ancova", ".csv", date=TRUE))

###############################################
#### Get slopes and intercepts for plotting ###
###############################################

# Note that this calculates interaction effects for all variables
# (even though min temp isn't significant)

mod<- NULL
mod2 <- NULL
out <- NULL
slope.epi <- list()
inter.epi <- list()
slope.ter <- list()
inter.ter<- list()
r.squared <- list()
pval <- list ()
for (i in 1:length(resp_vars)) {
  mod <- lm(formula(paste(resp_vars[i], " ~ habit/el + 0", sep="")), data=daily_values)
  mod2 <- lm(formula(paste(resp_vars[i], " ~ habit*el", sep="")), data=daily_values)
  out <- coef(mod)
  slope.epi[[i]] <- out[3]
  inter.epi[[i]] <- out[1]
  slope.ter[[i]] <- out[4]
  inter.ter[[i]] <- out[2]
  r.squared[[i]] <- summary(mod2)$r.squared
  pval[[i]] <- summary(mod2)$coefficients[4,4]
}

slopes <- data.frame(r.squared = unlist(r.squared), slope.epi = unlist(slope.epi), inter.epi=unlist(inter.epi), slope.ter=unlist(slope.ter), inter.ter=unlist(inter.ter), pval = unlist(pval))
rownames(slopes) <- resp_vars

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
partA <- make_scatter_plot (grand_means, slopes, "el", "mean_temp", "Mean Temp. (C)")
partA <- partA + blank_x.theme

partB <- make_scatter_plot (grand_means, slopes, "el", "sd_temp", "SD Temp. (C)")
partB <- partB + blank_x.theme

partC <- make_scatter_plot (grand_means, slopes, "el", "min_RH", "Minimum RH (%)")

partD <- make_scatter_plot (grand_means, slopes, "el", "sd_RH", "SD RH (%)")

# assemble the figure parts together
climate_plot <- plot_grid(partA, partB, partC, partD, nrow=2, ncol=2, labels = NULL, align = "hv")

# write out pdf
save_plot("microclimate.pdf", climate_plot, base_height = 4, base_width = 6)
