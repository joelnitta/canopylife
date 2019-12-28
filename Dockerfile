FROM rocker/geospatial:3.6.0

ARG DEBIAN_FRONTEND=noninteractive

# Install R packages with renv

COPY ./renv.lock ./

COPY ./renv_restore.R ./

RUN mkdir renv

RUN Rscript renv_restore.R

# Modify Rprofile.site so R loads renv library by default

RUN echo '.libPaths("/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Install latex packages with tinytex

COPY ./install_latex.R ./

RUN Rscript install_latex.R

# Set working dir

WORKDIR /home/rstudio/
