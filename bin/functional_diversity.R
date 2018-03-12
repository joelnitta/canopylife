# functional_diversity

# Run multivariate and univariate analyses of functional and phylogenetic diversity
# of epiphytic and terrestrial ferns from Moorea

# load packages
library(FD)
library(picante)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(mooreaferns)

# set working directory
setwd(here::here())

# define helper functions -------------------------------------------------

# wrapper to calculate mpd, mntd for epiphytic vs terrestrial communities
comstruc.habit <- function (comm, phy, traits) {
  
  # check for mismatches/missing species between phy and comm
  combined <- match.phylo.comm(phy, comm)
  phy <- combined$phy
  comm <- combined$comm
  
  ### split communities into epiphytic / terrestrial taxa, but use whole fern community as null model
  # get list of epiphytic and terrestrial taxa
  epi_species <- rownames(traits)[traits$habit == "epiphytic"]
  ter_species <- rownames(traits)[traits$habit == "terrestrial"]
  
  # define exclusively epiphytic or terrestrial communities 
  comm.epi <- comm
  comm.epi[,colnames(comm.epi) %in% ter_species] <- 0
  rownames(comm.epi) <- paste(rownames(comm.epi), "E", sep="_")
  
  comm.ter <- comm
  comm.ter[,colnames(comm.ter) %in% epi_species] <- 0
  rownames(comm.ter) <- paste(rownames(comm.ter), "T", sep="_")

  comm <- rbind(comm.epi, comm.ter)
  
  ### run community structure analysis ###
  
  # Look at phylo structure of each plot individually, compared to null (standard effect size). Null model shuffles tips of phylogeny in plots, keeps freq and richness the same
  # Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) indicate phylogenetic evenness.
  # Negative SES values and low quantiles (mpd.obs.p < 0.05) indicate phylogenetic clustering.
  # SES values of 0 are for trees with species spread randomly across the tree.
  # mpd.obs.z is the standardized effect size of mpd vs. null communities (equivalent to -NRI)
  
  # run ses.mpd
  mpd.out <- ses.mpd(comm, cophenetic(phy), null.model = "phylogeny.pool" , abundance.weighted = TRUE, runs = 999, iterations = 1000)
  # run ses.mntd
  mntd.out <- ses.mntd(comm, cophenetic(phy), null.model = "phylogeny.pool" , abundance.weighted = TRUE, runs = 999, iterations = 1000)
  
  # merge results
  mpd.out <-  mpd.out[, c("ntaxa", "mpd.obs", "mpd.obs.z", "mpd.obs.p")]
  mntd.out <-  mntd.out[, c("mntd.obs", "mntd.obs.z", "mntd.obs.p")]
  comstruc.out <- merge(mpd.out,mntd.out,by="row.names")
  colnames(comstruc.out) <- gsub("Row.names", "site", colnames(comstruc.out))
  comstruc.out$habit <- ""
  comstruc.out$habit[grep("_E", comstruc.out$site)] <- "epiphytic"
  comstruc.out$habit[grep("_T", comstruc.out$site)] <- "terrestrial"
  comstruc.out$site <- gsub("_E", "", comstruc.out$site)
  comstruc.out$site <- gsub("_T", "", comstruc.out$site)
  comstruc.out$site <- as.character(comstruc.out$site)
  
  return (comstruc.out)
}

# wrapper function to run FD analyses on community trait data 
run_FD <- function (traits, comm, traits.use, habit, abun.weighted = TRUE, print_CWM = FALSE) {
  
  # trim traits
  traits <- traits[,colnames(traits) %in% traits.use ]
  
  # subset to only epiphytic or terrestrial
  if (habit == "terrestrial") {
    traits <- traits[traits$habit == "terrestrial",]
  } else if (habit == "epiphytic") {
    traits <- traits[traits$habit == "epiphytic",]
  } else {
    stop ("Must choose epiphytic or terrestrial for habit")
  }
  
  # drop habit
  traits$habit <- NULL
  
  # drop species with NAs for more than half of traits
  trait_cutoff <- ncol(traits)*0.5
  traits <- traits[rowSums(is.na(traits)) < trait_cutoff, ]
  
  # only include species with trait data and with at least one occurrence in all plots
  comm <- comm[, colnames(comm) %in% rownames(traits)]
  comm <- comm[complete.cases(comm), ]
  comm <- comm[, colSums(comm) > 0 ]
  
  # trim traits to only species occurring in plots
  traits <- traits[colnames(comm),]
  
  # check order of community data and metadata
  if (!(all.equal(colnames(comm), rownames(traits)))) {
    stop ("community and metadata don't match")
  }
  
  # run FD analysis (normal)
  out <- dbFD(x=traits, a=comm, corr = "lingoes", m="max", w.abun=abun.weighted, calc.CWM=print_CWM, stand.FRic = FALSE)
  
  out <- as.data.frame(out)
  
  # add columns for elevation, growth habit
  out$habit <- habit
  
  # add site column
  out$site <- rownames(out)
  rownames(out) <- NULL
  
  colnames(out) <- gsub("nbsp", "ntaxa_with_traits", colnames(out))
  
  # either include community weighted means or not based on user input
  if (print_CWM == TRUE) {
    cwm_names <- paste("CWM", colnames(traits), sep=".")
    out <- out[,c("site", "ntaxa_with_traits", "FRic", "FEve", "FDiv", "habit", cwm_names)]
  } else if (print_CWM == FALSE) {
    out <- out[,c("site", "ntaxa_with_traits", "FRic", "FEve", "FDiv", "habit")]
  }
  return (out)
}

# load data ---------------------------------------------------------------

### sites
# use only sites on Moorea
sites <- mooreaferns::sites[grep("Aorai", mooreaferns::sites$site, invert=TRUE), ]

### traits
# load untransformed traits
traits <- mooreaferns::fern_traits

# add rownames to traits, drop species
rownames(traits) <- traits$species
traits$species <- NULL

# define sporo and gameto traits
sporo_traits <-c("stipe", "length", "width", "rhizome", "dissection", "pinna", "sla")
gameto_traits <- c("morphotype", "glands", "hairs", "gemmae")

# for functional diversity metrics, only use sporo traits
traits.use <- c(sporo_traits, "habit")

### community data
comm <- mooreaferns::sporocomm

# trim community data to only moorea sites
comm <- comm[rownames(comm) %in% sites$site, ]

### tree
phy <- mooreaferns::fern_tree

# run FD analysis on untransformed community weighted means ---------------

# run functional diversity analyses
fd.ter <- run_FD(habit = "terrestrial", traits = traits, traits.use = traits.use, comm=comm, print_CWM = TRUE)
fd.epi <- run_FD(habit = "epiphytic", traits = traits, traits.use = traits.use, comm=comm, print_CWM = TRUE)

# make into one dataframe
compare.df <- rbind(fd.ter, fd.epi)

# trim "CWM." from column names
colnames (compare.df) <- gsub("CWM.", "", colnames(compare.df))

# trim columns
cwm <- compare.df [, colnames(compare.df) %in% c("site", "habit", sporo_traits)]

# run FD analysis on multivariate data ------------------------------------

# transform and scale traits
traits.trans <- transform_traits(traits)

# run functional diversity analyses
fd.ter <- run_FD(habit = "terrestrial", traits = traits.trans, traits.use = traits.use, comm=comm, print_CWM = FALSE)
fd.epi <- run_FD(habit = "epiphytic", traits = traits.trans, traits.use = traits.use, comm=comm, print_CWM = FALSE)

# make into one dataframe
func_div <- rbind(fd.ter, fd.epi)

rm (fd.epi, fd.ter, compare.df)


# run phylogenetic diversity analysis -------------------------------------

phy_div <- comstruc.habit (comm, phy, traits)

rm(comm)

# calculate mean microclimate values --------------------------------------

source("bin/microclimate_means.R")

# merge results -----------------------------------------------------------

# compare.df includes all sites
compare.df <- left_join(phy_div, func_div) %>% left_join(cwm) %>% left_join(sites)

# compare.env includes ONLY sites with env. data
compare.env <- left_join(grand_means, compare.df)

# convert habit to factor
compare.df$habit <- factor(compare.df$habit, levels =c("terrestrial", "epiphytic"))
compare.env$habit <- factor(compare.env$habit, levels =c("terrestrial", "epiphytic"))

rm (cwm, func_div, phy_div, grand_means, daily_values)

# function to conduct t-test ----------------------------------------------

# use a two-sided test to check for differences between the groups
run_t_test <- function (input_data, resp_var) {
  ter_vals <- input_data[input_data$habit == "terrestrial", resp_var]
  epi_vals <- input_data[input_data$habit == "epiphytic", resp_var]
  out <- t.test(x = ter_vals, y = epi_vals, alternative = "two.sided", na.action="na.omit")$p.value
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

# functions to format data for plotting -----------------------------------

# function to run linear model, fill in plot and text data frames with model observations, parameters
fill_df <- function (habit, response, indep_var, compare.data, model_name, plot_data, text_data) {

  # choose best model based on subsetted data
  model <- lm (formula=paste(response, indep_var, sep=" ~ "), data = compare.data[compare.data$habit == habit,])
  # fill in text data frame
  text_data$r.squared <- summary(model)$r.squared
  text_data$p <- as.numeric(summary(model)$coefficients[,4][2])
  text_data$habit <- habit
  text_data$variable <- model_name
  
  # make vector of model points. start with vector of NA values, then fill in appropriately
  model_vector <- as.numeric(rep(NA, nrow(plot_data)))
  names(model_vector) <- rownames(plot_data)
  model_vector[names(fitted(model))] <- fitted(model)
  
  # add the model points vector to the plot data df (as the last column)
  plot_data$fitted <- model_vector
  # rename the last column with the model name
  colnames(plot_data)[length(colnames(plot_data))] <- model_name
  
  results <- list(plot_data, text_data)
  return(results)
}

# wrap the fill_df function in another function to automate it
# input_compare.df is the df with results of FD analysis
# y is vector of response variables (either functional diveristy or community weighted means)
# x is elevation or min RH
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
  
  # format r sq column for plotting with italics (need to use parse=TRUE in ggplot)
  # also save original r sq value for MS
  text.data$r.squared.plain <- text.data$r.squared
  text.data$r.squared <- paste("italic(R)^2 ==", round(text.data$r.squared, 2), sep=" " )
  
  output <- list(plot.data, text.data)
  return (output)
}

# function to calculate mean, sd of resp vars for barplot, add asterisks from significant t-tests
make_barplot_means <- function (input_data, input_resp_vars, t_test_results) {
  output <- list()
  for (i in 1:length(input_resp_vars)) {
    # run summary SE to calculate mean, sd, se. output is dataframe.
    output[[i]] <- Rmisc::summarySE(data = input_data, measurevar = input_resp_vars[i], groupvars = "habit", na.rm=TRUE)
    # save the name of the resp var
    output[[i]]$resp_var <- input_resp_vars[i]
    # need to make sure each df in list has same column names, or rbind.fill won't work
    colnames(output[[i]]) <- c("habit", "N", "mean", "sd", "se", "ci", "resp_var")
  }
  # for loop generates a list of dataframes. use rbind.fill to collapse them into one.
  output <- plyr::rbind.fill(output)
  
  #### add asterisks for t-test results
  # first convert to dataframe for joining with oupput
  t_test_results <- data.frame(resp_var = names(t_test_results), pval = t_test_results)
  # join the two dataframes. p values will be duplicated between epi and ter for each resp var, but that's OK
  output <- plyr::join (output, t_test_results, by = "resp_var", match ="first")
  # delete pvalues for epis (only need one pvalue for each comparison)
  output$pval[output$habit == "epiphytic"] <- ""
  output$pval <- as.numeric(output$pval)
  
  # make asterisks based on pval bins
  output$asterisk <- ""
  output$asterisk[output$pval < 0.05] <- "*"
  output$asterisk[output$pval < 0.01] <- "**"
  output$asterisk[output$pval < 0.001] <- "***"
  
  return (output)
}

# plot functions ----------------------------------------------------------

# plot function that adds lines based on significance (with legend removed)
make_scatter_plot <- function (input.plot.data, input.text.data, indep_var, dep_var, xlabel, ylabel) {
  p <- ggplot (data = input.plot.data, aes_string(x = indep_var, y = dep_var, color = "habit")) +
    geom_point(size = 2) +
    standard_theme +
    theme(legend.position = "none") +
    scale_colour_manual(values = cols) +
    xlab(xlabel) +
    ylab(ylabel)
  
  # get actual y and x axis ranges
  yrange <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
  xrange <- ggplot_build(p)$layout$panel_ranges[[1]]$x.range
  
  # add text in upper right hand corner for R2 value, with one lower than the other if there are two values
  if ( (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "terrestrial"] < 0.05) &
       (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "epiphytic"] < 0.05) ) {
    
    p <- p +
      geom_line(aes_string(x = indep_var, y = paste("ter", dep_var, "model", sep="."))) +
      geom_text(x=0.99*xrange[2],
                y=0.99*yrange[2],
                aes(label=r.squared),
                data=subset(input.text.data, habit=="terrestrial" & variable==dep_var),
                parse=TRUE,
                hjust=1.0,
                vjust=1.0,
                show.legend=FALSE)

    p <- p + geom_line(aes_string(x = indep_var, y = paste("epi", dep_var, "model", sep="."))) +
      geom_text(x=0.99*xrange[2],
                y=0.88*yrange[2],
                aes(label=r.squared), 
                data=subset(input.text.data, habit=="epiphytic" & variable==dep_var),
                parse=TRUE,
                hjust=1.0,
                vjust=1.0,
                show.legend=FALSE)
      
  } else if ( (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "terrestrial"] < 0.05) &
              (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "epiphytic"] > 0.05)) {
    p <- p +
      geom_line(aes_string(x = indep_var, y = paste("ter", dep_var, "model", sep="."))) +
      geom_text(x=0.99*xrange[2],
                y=0.99*yrange[2],
                aes(label=r.squared),
                data=subset(input.text.data, habit=="terrestrial" & variable==dep_var),
                parse=TRUE,
                hjust=1.0,
                vjust=1.0,
                show.legend=FALSE)

  } else if ( (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "epiphytic"] < 0.05) &
              (input.text.data$p[input.text.data$variable == dep_var & input.text.data$habit == "terrestrial"] > 0.05)) {
    p <- p +
      geom_line(aes_string(x = indep_var, y = paste("epi", dep_var, "model", sep="."))) +
      geom_text(x=0.99*xrange[2],
                y=0.99*yrange[2],
                aes(label=r.squared),
                data=subset(input.text.data, habit=="epiphytic" & variable==dep_var),
                parse=TRUE,
                hjust=1.0,
                vjust=1.0,
                show.legend=FALSE)
  }
  
  return(p)
}

# plot total means as barplot
make_bar_plot <- function (input.plot.data, selected_resp_var) {
  
  # use only the data for the selected response variable
  subsetted.data <- input.plot.data[input.plot.data$resp_var == selected_resp_var, ] 
  
  # extract the asterisk for the pvalue
  # should either be one asterisk and one one blank, or both blank. keep only single asterisk value or blank.
  pval <- subsetted.data$asterisk
  pval <- ifelse(any(grepl("\\*", pval)), pval[grepl("\\*", pval)], "")
  
  p <- ggplot (data = subsetted.data,
               aes(x = habit, y = mean, fill = habit)) +
    geom_bar (position=position_dodge(), stat="identity", width=0.75) +
    geom_errorbar (aes(ymax=mean+se, ymin=mean-se),
                   width = .2,
                   position = position_dodge(0.9)) +
    geom_signif(comparisons = list(c("terrestrial", "epiphytic")),
                annotations=pval,
                textsize=5,
                size=NA,
                vjust=1) +
    # optionally plot number of species with that trait
    # geom_text (aes(label = N), y=0, vjust=0, size=4) + 
    xlab("Growth Habit") +
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
    geom_signif(comparisons = list(c("terrestrial", "epiphytic")),
                map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05, " "=2),
                test="t.test",
                test.args=c(alternative="two.sided"),
                textsize=5,
                size=NA,
                vjust=1) +
    xlab("Growth Habit") +
    standard_theme +
    theme(legend.position = "none") +
    scale_fill_manual(values = cols)
  return (p)
}

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
# outputs list: 1 is plots, 2 is dataframe of linear model results
plot_scatter_dataset <- function (resp_vars, dataset, indep_var, xlabel, ylabels) {
  
  # generate dataframes for ggplot
  plot_data <- run_fill_df(dataset, resp_vars, indep_var)[[1]]
  text_data <- run_fill_df(dataset, resp_vars, indep_var)[[2]]
  
  # loop through plot function to make list of individual plots
  plots <- list()
  for (i in 1:length(resp_vars)) {
    plots[[i]] <- make_scatter_plot (plot_data, text_data, indep_var, resp_vars[i], xlabel, ylabels[i])
  }
  
  results <- list(plots, text_data, plot_data)

  return(results)
  
}

# plotting: multivariate analysis ---------------------------------------------

#### general plot settings

# choose response variables for overall plot (here, functional diversity metrics)
resp_vars <- c("ntaxa", "mpd.obs.z", "mntd.obs.z", "FRic", "FEve", "FDiv")

# make list of y axis labels
my.ylabels <- c("Richness", "MPD", "MNTD", "FRic", "FEve", "FDiv")

# set colors: either default ggplot or green/brown
# cols <- gg_color_hue(2) # red/ter is first color, blue/epi is second
cols <- brewer.pal(9, "Set1")[c(7,3)] #green / brown from color brewer

#### make boxplots

# make dataframe for t-test values
# run t-test on entire dataset
t.test.out <- run_t_test_dataset (compare.df, resp_vars)

# make df of barplot data (don't need the means etc, but use the asterisks for the boxplot)
multvar.t.results <- make_barplot_means (compare.df, resp_vars, t.test.out)

# generate boxplots as list
boxplots <- plot_box_dataset (compare.df, resp_vars, multvar.t.results)

# modifications for boxplots
# strip y axis titles and scale, since these are same across all plots in a row
 for (i in 1:length(boxplots)) {
   boxplots[[i]] <- boxplots[[i]] +
     theme(axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.text.y = element_blank())
 }

# strip x axis from all but last
for (i in 1:(length(boxplots)-1)) {
  boxplots[[i]] <- boxplots[[i]] +
    blank_x.theme
}

#### generate scatter plots along gradient
# left side scatter plots have x-axis = min RH, right side scatter plots have x-axis = elevation
left_plots <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)[[1]]
right_plots <- plot_scatter_dataset (resp_vars, compare.df, "el", xlabel="Elevation (m)", ylabels=my.ylabels)[[1]]

# also output linear model results in dataframe for MS
multvar.env.results <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)[[2]]
multvar.el.results <- plot_scatter_dataset (resp_vars, compare.df, "el", xlabel="Elevation (m)", ylabels=my.ylabels)[[2]]

### modifications for scatter plots

# run out of space with automatic x axis for elevation, so set manually
for (i in 1:(length(right_plots))) {
  right_plots[[i]] <- right_plots[[i]] +
    scale_x_continuous(breaks=c(0,250,500,750,1000))
}

# the left side plots (by min. humidity) will have a smaller range of values than the right side plots (by elevation)
# to make the y axis match, need to expand range of each plot to match right side values
for (i in 1:(length(left_plots))) {
  left_plots[[i]] <- left_plots[[i]] +
    expand_limits(y = c(
      min(c(layer_scales(left_plots[[i]])$y$range$range[1], layer_scales(right_plots[[i]])$y$range$range[1], layer_scales(boxplots[[i]])$y$range$range[1])),
      max(c(layer_scales(left_plots[[i]])$y$range$range[2], layer_scales(right_plots[[i]])$y$range$range[2], layer_scales(boxplots[[i]])$y$range$range[2])) )
    )
}

for (i in 1:(length(right_plots))) {
  right_plots[[i]] <- right_plots[[i]] +
    expand_limits(y = c(
      min(c(layer_scales(left_plots[[i]])$y$range$range[1], layer_scales(right_plots[[i]])$y$range$range[1], layer_scales(boxplots[[i]])$y$range$range[1])),
      max(c(layer_scales(left_plots[[i]])$y$range$range[2], layer_scales(right_plots[[i]])$y$range$range[2], layer_scales(boxplots[[i]])$y$range$range[2])) )
    )
}

for (i in 1:(length(boxplots))) {
  boxplots[[i]] <- boxplots[[i]] +
    expand_limits(y = c(
      min(c(layer_scales(left_plots[[i]])$y$range$range[1], layer_scales(right_plots[[i]])$y$range$range[1], layer_scales(boxplots[[i]])$y$range$range[1])),
      max(c(layer_scales(left_plots[[i]])$y$range$range[2], layer_scales(right_plots[[i]])$y$range$range[2], layer_scales(boxplots[[i]])$y$range$range[2])) )
    )
}

# modifications for scatterplots
# strip y axis titles, ticks, and text for right-hand side of part A
for (i in 1:length(right_plots)) {
  right_plots[[i]] <- right_plots[[i]] +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

# remove x axis labels from all plots but the two on the bottom
for (i in 1:(length(left_plots)-1)) {
  left_plots[[i]] <- left_plots[[i]] +
    blank_x.theme
  right_plots[[i]] <- right_plots[[i]] +
    blank_x.theme
}

### use cowplot to assemble parts into final figure

# this gets all subplots aligned correctly, but results in too much whitespace
# no simple coding solution available, so trim this out with a graphics editor after
plots <- c(left_plots, right_plots, boxplots)

# plot grid will go left to right, top to bottom, so rearrange list of plots to match
resp_vars_sort <- c("1_ntaxa", "2_mpd.obs.z", "3_mntd.obs.z", "4_FRic", "5_FEve", "6_FDiv")
plot_list <- data.frame(type=c(rep("a_rh", 6), rep("b_el", 6), rep("c_box", 6)), resp_var = rep(resp_vars_sort, 3))
plot_list <- plot_list[order(plot_list$resp_var, plot_list$type), ]
plots <- plots[as.numeric(rownames(plot_list))]

fdiv_plots <- plot_grid(plotlist=plots, labels=NULL, align = "hv", ncol=3, nrow=6)
save_plot(fdiv_plots, file = "fdiv_plots.pdf", ncol=3, nrow=6, base_height = 2, base_width=2.5)

# plotting: univariate analysis -----------------------------------------------

### overall plot settings
# choose response variables
resp_vars <- c("stipe", "width", "rhizome", "dissection", "pinna", "sla")

# y axis labels
my.ylabels <- c("Stipe L. (cm)", "Frond W. (cm)", "Rhizome Dia. (cm)", "Dissection", "No. Pinna Pairs", "SLA (m^2 / kg)")

# set colors: green / brown from color brewer
cols <- brewer.pal(9, "Set1")[c(7,3)]

### generate plot dataframes, make plots

# for barplots, calculate overall means for untransformed, selected response variables (NOT community values)
traits.test <- traits[c("habit", resp_vars)]

# first, run t-test on all response variables
t.test.out <- run_t_test_dataset (traits.test, resp_vars)

# make barplot (have to generate dataframe first)
univar.t.results <- make_barplot_means (traits.test, resp_vars, t.test.out)

### make barplots
barplots <- plot_bar_dataset(univar.t.results, resp_vars)

# modifications for barplots
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

### make scatter plots
# left side scatter plots have x-axis = min RH, right side scatter plots have x-axis = elevation
left_plots <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)[[1]]
right_plots <- plot_scatter_dataset (resp_vars, compare.df, "el", xlabel="Elevation (m)", ylabels=my.ylabels)[[1]]

# also output linear model results in dataframe for MS
univar.env.results <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)[[2]]
univar.el.results <- plot_scatter_dataset (resp_vars, compare.env, "min_RH", xlabel="Minimum RH (%)", ylabels=my.ylabels)[[2]]

### make various adjustments to plots

# reformat SLA plot y-axis title to use superscript
left_plots[[which(resp_vars=="sla")]] <-
  left_plots[[which(resp_vars=="sla")]] +
  labs(y = expression(SLA~(m^2 / kg)))

# run out of space with automatic x axis for elevation, so set manually
for (i in 1:(length(right_plots))) {
  right_plots[[i]] <- right_plots[[i]] +
    scale_x_continuous(breaks=c(0,250,500,750,1000))
}

# remove x axis labels from all plots but the two on the bottom
# plot lists are same length, so can combine
for (i in 1:(length(left_plots)-1)) {
  left_plots[[i]] <- left_plots[[i]] +
    blank_x.theme
  right_plots[[i]] <- right_plots[[i]] +
    blank_x.theme
}

# the left side plots (by min. humidity) will have a smaller range of values than the right side plots (by elevation)
# to make the y axis match, need to expand range of each plot to match right side values
for (i in 1:(length(left_plots))) {
  left_plots[[i]] <- left_plots[[i]] +
    expand_limits(y = c(
      layer_scales(right_plots[[i]])$y$range$range[1],
      layer_scales(right_plots[[i]])$y$range$range[2])
    )
}

# modifications for scatterplots
# strip y axis titles, ticks, and text for right-hand scatter plots
for (i in 1:length(right_plots)) {
  right_plots[[i]] <- right_plots[[i]] +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

### use cowplot to merge all the subplots together, with equal size for each plot area

# this gets all subplots aligned correctly, but results in too much whitespace
# no simple coding solution available, so trim out whitespace with a graphics editor afterwards
plots <- c(left_plots, right_plots, barplots)

# plot grid will go left to right, top to bottom, so rearrange list of plots to match
resp_vars_sort <- paste(1:length(resp_vars), resp_vars, sep="_")
plot_list <- data.frame(type=c(rep("a_rh", 6), rep("b_el", 6), rep("c_bar", 6)), resp_var = rep(resp_vars_sort, 3))
plot_list <- plot_list[order(plot_list$resp_var, plot_list$type), ]
plots <- plots[as.numeric(rownames(plot_list))]

cwm_plots <- plot_grid(plotlist=plots, labels=NULL, align = "hv", ncol=3, nrow=6)
save_plot(cwm_plots, file = "cwm_plots.pdf", ncol=3, nrow=6, base_height = 2, base_width=2.5)
