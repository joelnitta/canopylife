#! /bin/bash
# Take a snapshot then build the image
docker run --rm -e DISABLE_AUTH=true -v /Users/joelnitta/R/canopylife:/home/rstudio/project rocker/geospatial:3.6.0 bash /home/rstudio/project/install_packages.sh
DOCKER_BUILDKIT=1 docker build --progress=plain --no-cache --secret id=pat,src=.secret . -t joelnitta/canopylife:3.6.0
rm .Rprofile
