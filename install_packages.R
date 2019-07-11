# This script is meant to be called by `build_image.sh` to
# install packages to the project library and make packrat/packrat.lock
# It should be run from the rocker/verse:3.6.0 docker image.
#
# This only is meant for the *initial* creation of packrat.lock;
# that file may be subsequently updated as packages are added/updated.
#
# For more info on installing R packages to docker images with
# packrat, see https://www.joelnitta.com/post/docker-and-packrat/

# Load secret GitHub PAT for installing private packages
Sys.setenv(GITHUB_PAT = readLines(".secret"))

### Initialize packrat ###

# Don't let packrat try to find
# packages to install itself.

install.packages("packrat", repos = "https://cran.rstudio.com/")
packrat::init(
  infer.dependencies = FALSE,
  enter = TRUE,
  restart = FALSE)

### Setup repositories ###

# Install packages that install packages.
install.packages("BiocManager", repos = "https://cran.rstudio.com/")
install.packages("remotes", repos = "https://cran.rstudio.com/")

# Specify repositories so they get included in
# packrat.lock file.
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
  "cowplot",
  "drake",
  "future",
  "gghighlight",
  "ggrepel",
  "ggridges",
  "glue",
  "here",
  "janitor",
  "kableExtra",
  "phytools",
  "picante",
  "scico",
  "tictoc",
  "tidyverse",
  "vegan",
  "viridis",
  "visNetwork"
  )

install.packages(cran_packages)

### Install bioconductor packages ###
bioc_packages <- c(
  "ggtree"
)

BiocManager::install(bioc_packages)

### Install github packages ###
# (Need to install dev version of treeio until Bioconductor
# catches up to 1.9.1 because of ape::root scoping bug.)
github_packages <- c(
  "GuangchuangYu/treeio",
  "joelnitta/jntools",
  "joelnitta/mooreaferns",
  "rstudio/gt",
  "thomasp85/patchwork"
)

remotes::install_github(github_packages)

### Take snapshot ###

packrat::snapshot(
  snapshot.sources = FALSE,
  ignore.stale = TRUE,
  infer.dependencies = FALSE)
