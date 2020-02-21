# Install packages to a docker image with renv
#
# This script should be run from the rocker/geospatial:3.6.0 docker image.
# This only is meant for the *initial* creation of renv.lock;
# that file may be subsequently updated as packages are added/updated.
#
# To update a single package after the image is made...
# Launch this container: `docker-compose up -d`
# Within container, start R: `R`
# Within the R session, specify repos:
# ```
# my_repos <- BiocManager::repositories()
# my_repos["CRAN"] <- "https://cran.rstudio.com/"
# options(repos = my_repos)
# ```
# Update(install) your package of choice: `install.packages("whatever")`
# Snapshot specifying the renv library: `renv::snapshot(type = "simple", library = "/renv")`

### Initialize renv ###

# Initialize renv, but don't let it try to find packages to install itself.
install.packages("remotes", repos = "https://cran.rstudio.com/")
# Use dev version with most recent bug fixes
remotes::install_github("rstudio/renv")

renv::consent(provided = TRUE)

renv::init(
  bare = TRUE,
  force = TRUE,
  restart = FALSE)

renv::activate()

### Setup repositories ###

# Install packages that install packages.
install.packages("BiocManager", repos = "https://cran.rstudio.com/")
# (Need to do this again because now we're in a fresh renv project)
install.packages("remotes", repos = "https://cran.rstudio.com/")

# Set repos.
my_repos <- BiocManager::repositories()
my_repos["CRAN"] <- "https://cran.rstudio.com/"
options(repos = my_repos)

### Install CRAN packages ###
cran_packages <- c(
  "FD",
  "FactoMineR",
  "RColorBrewer",
  "assertr",
  "assertthat",
  "bookdown",
  "broom",
  "caper",
  "checkr",
  "conflicted",
  "corrr",
  "cowplot",
  "drake",
  "future",
  "ggrepel",
  "ggridges",
  "ggtree",
  "glue",
  "here",
  "janitor",
  "kableExtra",
  "knitr",
  "mgcv",
  "phytools",
  "plantecophys",
  "picante",
  "pryr",
  "scico",
  "sme",
  "sp",
  "spaMM"
  "spdep",
  "taxize",
  "tictoc",
  "tidyverse",
  "vegan",
  "viridis",
  "visNetwork",
  "tinytex"
)

install.packages(cran_packages)

### Install github packages ###
github_packages <- c(
  "adletaw/captioner",
  "joelnitta/jntools",
  "rstudio/gt",
  "thomasp85/patchwork",
  "beckyfisher/FSSgam_package"
)

remotes::install_github(github_packages)

### Take snapshot ###

renv::snapshot(type = "simple")
