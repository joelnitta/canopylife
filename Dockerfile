# syntax = docker/dockerfile:1.0-experimental
FROM rocker/verse:3.6.0

ARG DEBIAN_FRONTEND=noninteractive

#######################################
### Install R packages with packrat ###
#######################################

# Only run this after making packrat/packrat.lock by
# running install_packages.R

COPY ./packrat/packrat.lock packrat/

COPY packrat_restore.R .

# RUN install2.r packrat
RUN Rscript -e 'install.packages("packrat", repos = "https://cran.rstudio.com/")'

RUN --mount=type=secret,id=pat \
Rscript packrat_restore.R `cat /run/secrets/pat`

# Modify Rprofile.site so R loads packrat library by default
RUN echo '.libPaths("/packrat/lib/x86_64-pc-linux-gnu/3.6.0")' >> /usr/local/lib/R/etc/Rprofile.site

WORKDIR /home/rstudio/
