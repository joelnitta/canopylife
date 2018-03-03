#!bin/bash

# deletes intermediate files. 
# these can't be recovered, so use with care!

# various intermediate latex files
rm -f *.aux
rm -f *.bbl
rm -f *.blg
rm -f *.dvi
rm -f *.log

# R output
rm -f results.RData
rm -f PCA.pdf
rm -f microclimate.pdf

# etc
rm -f readme.html
