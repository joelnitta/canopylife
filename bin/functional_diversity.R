library(vegan)
library(picante)
library(ade4)
library(FD)
library(ggplot2)
library(reshape)
library(xlsx)
library(FactoMineR)
library(gridExtra)
library(cowplot)
library(RColorBrewer)

library(gridExtra)
library(grid)
library(gtable)
library(plyr)

setwd("/Users/joelnitta/R/moorea/")
source("bin/my_funcs.R")
source("bin/comstruc.R")
source("bin/my_themes.R")
source("bin/get_hobo_summary.R")
source("bin/traits/get_traits.R")
source("bin/traits/run_FD.R")
source("bin/traits/comstruc_habit.R")
source("bin/microclimate/hobo_threshold.R")

# make plot of various community functional trait diversity metrics against elevation
# also run PCA on Moorea environmental data, plot against functional diveristy

# set options
# traits.select <- select.list(c("gameto", "sporo", "both"), title="Select traits")
traits.select <- "sporo"
# abundance <- select.list (c("abun_weighted", "not abun_weighted"), title="Select abundance weighting for FD analysis")
abundance <- "abun_weighted"

###############################
#### environmental analysis ###
###############################

# load site data
lat_long_el <- get.my.env("data/sites.csv", site = "moorea", both.gen = FALSE)

# load raw hobo data
hobo.raw <- read.csv("data/hobo_moorea_aorai_2-24-15.csv",header=T)

# threshold at chosen value for 2 hr
thresh <- 80
threshold.df <- calc_thresh(hobo.raw, thresh, 8, lat_long_el)

# run env PCA
source("bin/traits/chp3_env_PCA.R")

rm (env, env.for.pca, env.scaled, hobo.raw, variables.pca)

#################################################################
### run FD analysis on untransformed community weighted means ###
#################################################################

# load untransformed traits
traits <- get.my.traits(log.trans = FALSE, scale.traits = FALSE)

sporo_traits <-c("stipe", "length", "width", "rhizome", "dissection", "pinna", "SLA")
gameto_traits <- c("morphotype", "glands", "hairs", "gemmae")

traits_use <- c(sporo_traits, "habit")

fd.ter <- run_FD("ter", abun.setting = ifelse(abundance == "abun_weighted", TRUE, FALSE), traits, traits_use, print_CWM = TRUE)
fd.epi <- run_FD("epi", abun.setting = ifelse(abundance == "abun_weighted", TRUE, FALSE), traits, traits_use, print_CWM = TRUE)

# add site_height column (uniquely identify epi or ter and site name) to merge with datalogger PCA
fd.ter$site_height <- paste(rownames(fd.ter), "ter", sep = "_")
fd.epi$site_height <- paste(rownames(fd.epi), "epi", sep = "_")

# make into one dataframe
compare.df <- rbind(fd.ter, fd.epi)

# convert habit to factor
compare.df$habit <- factor(compare.df$habit, levels =c("ter", "epi"))
levels(compare.df$habit) <- c("terrestrial", "epiphytic")

# trim "CWM." from column names
colnames (compare.df) <- gsub("CWM.", "", colnames(compare.df))

# trim columns
cwm.allsites <- compare.df [, colnames(compare.df) %in% c("site_height", sporo_traits)]
# compare.cwm <- merge(compare.cwm, env.pca[,c("Dim.1", "site_height")], by = "site_height")
# compare.cwm <- merge(compare.cwm, threshold.df, by.x = "site_height", by.y = "site")

rm (compare.df, fd.epi, fd.ter, traits)

############################################
### run FD analysis on multivariate data ###
###########################################

# load transformed, scaled traits
traits <- get.my.traits(log.trans = TRUE, scale.traits = TRUE)

sporo_traits <-c("stipe", "length", "width", "rhizome", "dissection", "pinna", "SLA")
gameto_traits <- c("morphotype", "glands", "hairs", "gemmae")

traits_use <- c(sporo_traits, "habit")

fd.ter <- run_FD("ter", abun.setting = ifelse(abundance == "abun_weighted", TRUE, FALSE), traits, traits_use)
fd.epi <- run_FD("epi", abun.setting = ifelse(abundance == "abun_weighted", TRUE, FALSE), traits, traits_use)

# add site_height column (uniquely identify epi or ter and site name) to merge with datalogger PCA
fd.ter$site_height <- paste(rownames(fd.ter), "ter", sep = "_")
fd.epi$site_height <- paste(rownames(fd.epi), "epi", sep = "_")

# make into one dataframe
compare.df <- rbind(fd.ter, fd.epi)

# convert habit to factor
compare.df$habit <- factor(compare.df$habit, levels =c("ter", "epi"))
levels(compare.df$habit) <- c("terrestrial", "epiphytic")

fdiv.allsites <- compare.df

rm (fd.epi, fd.ter, traits, compare.df)

##############################
### phylogenetic diversity ###
##############################

# load transformed, scaled traits
traits <- get.my.traits(log.trans = TRUE, scale.traits = TRUE)

comm_file <- set_input_files()[[1]]
tree_file <- set_input_files()[[2]]

comm <- get.my.comm(comm_file, site = "moorea", generation = "sporo", dataset = "full", abundance = "raw")
phy <- get.my.tree(tree_file)

phy.allsites <- comstruc.habit (comm, phy, traits, lat_long_el)
phy.allsites$site_height <- rownames(phy.allsites)
phy.allsites$site_height <- gsub("_T", "_ter", phy.allsites$site_height)
phy.allsites$site_height <- gsub("_E", "_epi", phy.allsites$site_height)
phy.allsites <- phy.allsites[,c("mpd.obs.z", "mpd.obs.p", "mntd.obs.z", "mntd.obs.p", "site_height")]

rm(comm, traits)

#########################################
#### merge results together for plotting
#########################################

# compare.df includes all sites
compare.df <- merge(fdiv.allsites, phy.allsites, by = "site_height")
compare.df <- merge(compare.df, cwm.allsites, by = "site_height")

# compare.env includes ONLY sites with env. data
compare.env <- merge(compare.df, env.pca[,c("Dim.1", "site_height")], by = "site_height")
compare.env <- merge(compare.env, threshold.df, by.x = "site_height", by.y = "site")

# optional: add min/mean/max/sd climate values to compare.env
df_to_add <- env.means
df_to_add$site_height <- paste(df_to_add$site, df_to_add$height, sep = "_")
df_to_add <- df_to_add[,grepl("max|mean|min|sd|site_height", colnames(df_to_add))]
compare.env <- merge(compare.env, df_to_add, by ="site_height")

rm (cwm.allsites, fdiv.allsites, phy.allsites, threshold.df, env.pca, lat_long_el, df_to_add )

#####################################
### t-test between the two groups ###
#####################################
# set up as list for plotting
# resp_vars <- c("nbsp", "mpd.obs.z", "mntd.obs.z", "FRic", "FEve", "FDiv")

run_t_test <- function (input_data, resp_var) {
  ter_vals <- input_data[input_data$habit == "terrestrial", resp_var]
  epi_vals <- input_data[input_data$habit == "epiphytic", resp_var]
  out <- t.test(x = ter_vals, y = epi_vals, alternative = "greater")$p.value
  return (out)
}

run_t_test_dataset <- function (input_data, input_resp_vars) {
  output <- list()
  for (i in 1:length(input_resp_vars)) {
    output[i] <- run_t_test (input_data, input_resp_vars[i])
  }
  names(output) <- input_resp_vars
  output <- unlist(output)
  return(output)
}

################################################
### check for significance with linear model ###
################################################

# this is now automatically done in plot function, but can run here if want
# to examine model output

# elevation dataset
# summary(lm(SLA~ el, data=compare.df[compare.df$habit == "terrestrial",]))
# summary(lm(FRic ~ el, data=compare.df[compare.df$habit == "epiphytic",]))

# environmental dataset
# summary(lm(FRic ~ el, data=compare.env[compare.env$habit == "terrestrial",]))
# summary(lm(FRic ~ el, data=compare.env[compare.env$habit == "epiphytic",]))

###############################
### data setup functions
###############################

# function to run linear model, fill in plot and text data frames with model observations, parameters
fill_df <- function (habit, response, input_indep_var, input_compare.data, model_name, plot_data, text_data) {
  plot.data <- plot_data
  text.data <- text_data
  
  # rename input variables (same name inside and outside the function causes problems!)
  indep_var <- input_indep_var
  compare.data <- input_compare.data
  
  # choose best model based on subsetted data
  model <- lm (formula=paste(response, indep_var, sep=" ~ "), data = compare.data[compare.data$habit == habit,])
  # fill in text data frame
  text.data$r.squared <- summary(model)$r.squared
  text.data$p <- as.numeric(summary(model)$coefficients[,4][2])
  text.data$habit <- habit
  text.data$variable <- model_name
  
  # make vector of model points. start with vector of NA values, then fill in appropriately
  model_vector <- as.numeric(rep(NA, nrow(plot.data)))
  names(model_vector) <- rownames(plot.data)
  model_vector[names(fitted(model))] <- fitted(model)
  
  # add the model points vector to the plot data df (as the last column)
  plot.data$fitted <- model_vector
  # rename the last column with the model name
  colnames(plot.data)[length(colnames(plot.data))] <- model_name
  
  results <- list(plot.data, text.data)
  return(results)
}

# wrap the fill_df function in another function to automate it
# input_compare.df is the df with results of FD analysis
# y is vector of response variables (either functional diveristy or community weighted means)
# x is elevation, number of events, or PC1
run_fill_df <- function (df, y, x) {
  
  # first plot.data is the same as the input comparative df
  plot.data <- df
  text.data <- data.frame(indep_variable = x)
  
  # fill in plot data using loop over fill_df function
  for (i in 1:length(y)) {
    plot.data <- fill_df("terrestrial", y[i], x, df, paste("ter", y[i], "model", sep ="."), plot.data, text.data)[[1]]
    plot.data <- fill_df("epiphytic", y[i], x, df, paste("epi", y[i], "model", sep ="."), plot.data, text.data)[[1]]
  }
  
  # do first text data alone
  text.data <- data.frame(indep_variable = x)
  text.data <- fill_df ("terrestrial", y[1], x, df, y[1], plot.data, text.data) [[2]]
  # now loop
  for (i in 1:length(y)) {
    text.data <- rbind(text.data, fill_df ("terrestrial", y[i], x, df, y[i], plot.data, text.data) [[2]])
    text.data <- rbind(text.data, fill_df ("epiphytic", y[i], x, df, y[i], plot.data, text.data) [[2]])
  }
  # git rid of redundant values
  text.data <- unique(text.data)
  
  # format r sq column for plotting
  text.data$r.squared <- paste("r2 = ", round(text.data$r.squared, 2), sep="" )
  
  output <- list(plot.data, text.data)
  return (output)
}

# function to calculate mean, sd of resp vars for barplot, add asterisks from significant t-tests
make_barplot_means <- function (input_data, input_resp_vars, t_test_results) {
  output <- list()
  for (i in 1:length(input_resp_vars)) {
    # run summary SE to calculate mean, sd, se. output is dataframe.
    output[[i]] <- summarySE(data = input_data, measurevar = input_resp_vars[i], groupvars = "habit", na.rm=TRUE)
    # save the name of the resp var
    output[[i]]$resp_var <- input_resp_vars[i]
    # need to make sure each df in list has same column names, or rbind.fill won't work
    colnames(output[[i]]) <- c("habit", "N", "mean", "sd", "se", "ci", "resp_var")
  }
  # for loop generates a list of dataframes. use rbind.fill to collapse them into one.
  output <- rbind.fill(output)
  
  #### add asterisks for t-test results
  # first convert to dataframe for joining with oupput
  t_test_results <- data.frame(resp_var = names(t_test_results), pval = t_test_results)
  # join the two dataframes. p values will be duplicated between epi and ter for each resp var, but that's OK
  output <- join (output, t_test_results, by = "resp_var", match ="first")
  # delete pvalues for epis (since it's one-sided, will only print asterisk above terrestrial)
  output$pval[output$habit == "epiphytic"] <- ""
  output$pval <- as.numeric(output$pval)
  
  # make asterisks based on pval bins
  output$asterisk <- ""
  output$asterisk[output$pval < 0.05] <- "*"
  output$asterisk[output$pval < 0.01] <- "**"
  output$asterisk[output$pval < 0.001] <- "***"
  
  # convert pvalues to asterisks
  return (output)
}

#####################
### plot functions ###
#####################

# plot function that adds lines based on significance (with legend removed)
make_scatter_plot <- function (input.plot.data, input.text.data, indep_var, dep_var, xlabel, ylabel) {
  p <- ggplot (data = input.plot.data, aes_string(x = indep_var, y = dep_var, color = "habit")) +
    geom_point(size = 2) +
    standard_theme +
    theme(legend.position = "none") +
    scale_colour_manual(values = cols) +
    xlab(xlabel) +
    ylab(ylabel)
  
  if (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "terrestrial"] < 0.05) {
    p <- p +
      geom_line(aes_string(x = indep_var, y = paste("ter", dep_var, "model", sep="."))) +
      geom_text(x=Inf, 
                y=Inf, 
                aes(label=r.squared), 
                data=subset(input.text.data, habit=="terrestrial" & variable==dep_var), 
                parse=FALSE, 
                hjust=1.0,
                vjust=1.0,
                #  vjust=c(1,2),
                show.legend=FALSE)
  }
  if (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "epiphytic"] < 0.05) {
    p <- p +
      geom_line(aes_string(x = indep_var, y = paste("epi", dep_var, "model", sep="."))) +
      geom_text(x=Inf, 
                y=Inf, 
                aes(label=r.squared), 
                data=subset(input.text.data, habit=="epiphytic" & variable==dep_var), 
                parse=FALSE, 
                hjust=1.0,
                vjust=1.0,
                # vjust=c(1,2),
                show.legend=FALSE)
  }
  
  return(p)
}

# plot total means as barplot
make_bar_plot <- function (input.plot.data, selected_resp_var) {
  p <- ggplot (data = input.plot.data[input.plot.data$resp_var == selected_resp_var, ], 
               aes(x = habit, y = mean, fill = habit)) +
    geom_bar (position=position_dodge(), stat="identity", width=0.75) +
    geom_errorbar (aes(ymax=mean+se, ymin=mean-se),
                   width = .2,
                   position = position_dodge(0.9)) +
    geom_text (aes(label = asterisk), size = 6) +
    geom_text (aes(label = N), y=0, vjust=-1) +
    standard_theme +
    theme(legend.position = "none") +
    scale_fill_manual(values = cols)
  return (p)
}

# plot total means as boxplot
make_box_plot <- function (input.plot.data, selected_resp_var, asterisk_data) {
  p <- ggplot (data = input.plot.data, 
               aes_string(x = "habit", y = selected_resp_var, fill = "habit")) +
    geom_boxplot(width=0.75) +
    geom_text (data = asterisk_data[asterisk_data$resp_var == selected_resp_var,],
               aes(x=habit, y=mean, label=asterisk),
               size=6) +
    standard_theme +
    theme(legend.position = "none") +
    scale_fill_manual(values = cols)
  return (p)
}

#notrun: test <- make_box_plot (compare.df, "mntd.obs.z", barplot_data)

# wrapper function to apply boxplot across whole dataset
plot_box_dataset <- function (input_barplot_data, input_resp_vars, input_asterisk_data) {
  output <- list()
  for (i in 1:length(input_resp_vars)) {
    output[[i]] <- make_box_plot (input_barplot_data, input_resp_vars[i], input_asterisk_data)
  }
  return(output)
}

# wrapper function to apply bar plot across whole dataset
plot_bar_dataset <- function (input_barplot_data, input_resp_vars) {
  output <- list()
  for (i in 1:length(input_resp_vars)) {
    output[[i]] <- make_bar_plot (input_barplot_data, input_resp_vars[i])
  }
  return(output)
}

# wrapper function to apply scatterplot func to an entire dataset
plot_scatter_dataset <- function (the_resp_vars, the_dataset, the_indep_var, xlabel, ylabels) {
  # generate dataframes for ggplot
  the_plot_data <- run_fill_df(the_dataset, the_resp_vars, the_indep_var)[[1]]
  the_text_data <- run_fill_df(the_dataset, the_resp_vars, the_indep_var)[[2]]
  
  # loop through plot function to make list of individual plots
  the_plots <- list()
  for (i in 1:length(the_resp_vars)) {
    the_plots[[i]] <- make_scatter_plot (the_plot_data, the_text_data, the_indep_var, the_resp_vars[i], xlabel, ylabels[i])
  }
  
  # # add ticks back to last plot (but have to remove legend again)
  # the_plots[[length(the_plots)]] <- the_plots[[length(the_plots)]] +
  #   standard_theme +
  #   theme(legend.position = "none")
  
  return(the_plots)
}

################################### 
### plotting: multivar analysis ###
###################################

###### general plot settings

# choose response variables for overall plot (here, functional diversity metrics)
resp_vars <- c("nbsp", "mpd.obs.z", "mntd.obs.z", "FRic", "FEve", "FDiv")

# make list of y axis labels
my.ylabels <- c("Richness (No. spp.)", "MPD", "MNTD", "FRic", "FEve", "FDiv")

# set colors: either default ggplot or green/brown
# cols <- gg_color_hue(2) # red/ter is first color, blue/epi is second
cols <- brewer.pal(9, "Set1")[c(7,3)] #green / brown from color brewer

###### make boxplots (fig part b)

# make dataframe for t-test values
# run t-test on entire dataset
t.test.results <- run_t_test_dataset (compare.df, resp_vars)
# make df of barplot data (don't need the means etc, but use the asterisks for the boxplot)
barplot_data <- make_barplot_means (compare.df, resp_vars, t.test.results)

# generate boxplots as list
boxplots <- plot_box_dataset (compare.df, resp_vars, barplot_data)

# modifications for boxplots
# strip y axis titles (but keep scale, because different from part A)
 for (i in 1:length(boxplots)) {
   boxplots[[i]] <- boxplots[[i]] +
     theme(axis.title.y = element_blank())
 }

# strip x axis from all but last
for (i in 1:(length(boxplots)-1)) {
  boxplots[[i]] <- boxplots[[i]] +
    blank_x.theme
}

# combine into one object for final figure
boxplots <- plot_grid(boxplots[[1]], boxplots[[2]], boxplots[[3]], boxplots[[4]], boxplots[[5]], boxplots[[6]], labels = NULL, ncol = 1, align = 'hv')

##### generate scatter plots along gradient (part A)
left.plots <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)
right.plots <- plot_scatter_dataset (resp_vars, compare.df, "el", xlabel="Elevation (m)", ylabels=my.ylabels)

# modifications for scatter plots
# DOESN'T WORK: column width won't be the same.
# # remove y axis titles for middle column
# for (i in 1:length(right.plots)) {
#   right.plots[[i]] <- right.plots[[i]] +
#     theme(axis.title.y = element_blank())
# }

# remove x axis labels from all plots but the two on the bottom
# plot lists are same length, so can combine
for (i in 1:(length(left.plots)-1)) {
  left.plots[[i]] <- left.plots[[i]] +
    blank_x.theme
  right.plots[[i]] <- right.plots[[i]] +
    blank_x.theme
}

# use cowplot to plot together
# first, make three cowplot grids of each set
left.plots <- plot_grid(plotlist=left.plots, labels = NULL, ncol = 1, align = 'hv')
right.plots <- plot_grid(plotlist=right.plots, labels = NULL, ncol = 1, align = 'hv')

# then, combine, save to pdf
# our goal size for a full-page figure is 8.5 x 6.5 inches
# the plots have a pretty big margin, so can cut this down in illustrator
# export at 16 x 11 inches, then remove .4 inches between each row
# this will produce 14 x 11 plot, then by 75% to 8.5 x 6.5

# pdf(file = make_filename("FDiv_plots", ".pdf", date=TRUE), height = 16, width = 11)
# plot_grid(left.plots, right.plots, boxplots, labels = c("A", "", "B"), ncol = 3, align='hv')
# dev.off()

# realized it needs to be 6in wide, not 6.5
# this size works for trimming down to 6 x 8.5
pdf(file = make_filename("FDiv_plots", ".pdf", date=TRUE), height = 11, width = 6.25)
plot_grid(left.plots, right.plots, boxplots, labels = c("A", "", "B"), ncol = 3, align='hv')
dev.off()

#################################
### plotting: univar analysis ###
#################################

### overall plot settings
# choose response variables
# NOTE: going to trim out width and length, but leave in one of them for now so we have a total of
# six traits, which will produce individual plots the same size as for the FDiv plot (also 6 panels tall)
# resp_vars <- c("stipe", "length", "width", "rhizome", "dissection", "pinna", "SLA")
# resp_vars <- c("stipe", "rhizome", "dissection", "pinna", "SLA")
resp_vars <- c("stipe", "width", "rhizome", "dissection", "pinna", "SLA")

#y axis labels
# my.ylabels <- c("Stipe L. (cm)", "Rhizome Dia. (cm)", "Dissection", "No. Pinna Pairs", "SLA (m2 / kg)")
# my.ylabels <- c("Stipe L. (cm)", "Frond L. (cm)", "Frond W. (cm)", "Rhizome Dia. (cm)", "Dissection", "No. Pinna Pairs", "SLA (m2 / kg)")
my.ylabels <- c("Stipe L. (cm)", "Frond W. (cm)", "Rhizome Dia. (cm)", "Dissection", "No. Pinna Pairs", "SLA (m2 / kg)")

# set colors: either default ggplot or green/brown
# cols <- gg_color_hue(2) # red/ter is first color, blue/epi is second
cols <- brewer.pal(9, "Set1")[c(7,3)] #green / brown from color brewer

### generate plot dataframes, make plots
# for barplots, calculate overall means for selected response variables (NOT community values)
# load untransformed traits, convert habit to factor
traits <- get.my.traits(log.trans = FALSE, scale.traits = FALSE)
traits <- traits[c("habit", resp_vars)]
traits$habit <- as.factor(traits$habit)
traits$habit <- mapvalues(traits$habit, from = c("0", "1"), to = c("terrestrial", "epiphytic"))

# first, run t-test on all response variables
t.test.results <- run_t_test_dataset (traits, resp_vars)

# make barplot (have to generate dataframe first)
barplot.data <- make_barplot_means (traits, resp_vars, t.test.results)

##### make barplots (figure part B)
barplots <- plot_bar_dataset(barplot.data, resp_vars)

# modifications for boxplots
# strip y axis titles (but keep scale, because different from part A)
for (i in 1:length(barplots)) {
  barplots[[i]] <- barplots[[i]] +
    theme(axis.title.y = element_blank())
}
# strip x axis from all but last
for (i in 1:(length(barplots)-1)) {
  barplots[[i]] <- barplots[[i]] +
    blank_x.theme
}

##### make scatter plots using environmental dataset in one step (figure part B)
left.plots <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)
right.plots <- plot_scatter_dataset (resp_vars, compare.df, "el", xlabel="Elevation (m)", ylabels=my.ylabels)

# modifications for scatterplots
# strip y axis titles (but keep scale, because different from left side)
# for (i in 1:length(right.plots)) {
#   right.plots[[i]] <- right.plots[[i]] +
#     theme(axis.title.y = element_blank())
# }

# remove x axis labels from all plots but the two on the bottom
# plot lists are same length, so can combine
for (i in 1:(length(left.plots)-1)) {
  left.plots[[i]] <- left.plots[[i]] +
    blank_x.theme
  right.plots[[i]] <- right.plots[[i]] +
    blank_x.theme
}

# use cowplot to plot together
# first, make three cowplot grids of each set
barplots <- plot_grid(plotlist = barplots, labels = NULL, ncol = 1, align = 'hv')
left.plots <- plot_grid(plotlist = left.plots, labels = NULL, ncol = 1, align = 'hv')
right.plots <- plot_grid(plotlist = right.plots, labels = NULL, ncol = 1, align = 'hv')

# then, combine
# pdf(file = make_filename("CWM_plots", ".pdf", date=TRUE), height = 16, width = 11)
# plot_grid(left.plots, right.plots, barplots, labels = c("A", "", "B"), ncol = 3)
# dev.off()

# then, combine
pdf(file = make_filename("CWM_plots", ".pdf", date=TRUE), height = 11, width = 6.25)
plot_grid(left.plots, right.plots, barplots, labels = c("A", "", "B"), ncol = 3)
dev.off()