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

# Load functions and plans.
source("code/functions.R")
source("code/plan.R")

# Make plan.
make(plan)
