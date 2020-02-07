# make.R 
#
# Master script for running analyses.
#
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Set working directory to project root.
setwd(here::here())

# Load packages.
source("code/setup.R")

# Update drake settings.
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Specify non-global environment
# to get around taxize and captioner modifying global env
# (cf https://github.com/ropensci/drake/issues/749)
envir <- new.env(parent = globalenv())

# Load functions and plans.
source("code/functions.R", local = envir)
source("code/plan.R")

# Load cache
canopy_cache <- new_cache("canopy_cache")

# Set seed for reproducibility.
set.seed(1954)

# Make plan.
make(plan, cache = canopy_cache, envir = envir)
