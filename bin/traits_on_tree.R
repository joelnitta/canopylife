# traits_on_tree

# Plot sporophyte and gametophyte traits on phylogenetic tree

# load packages
library(ape)
library(scales)
library(RColorBrewer)
library(mooreaferns)

# axisGeo from phyloch is used to make geological time scale
# phyloch must be installed from source from the developer's website
# http://www.christophheibl.de/Rpackages.html
# download "package source", unzip, then run: install.packages("path_to_package", repos = NULL, type="source")
library(phyloch)

# set working directory
setwd(here::here())

# helper functions --------------------------------------------------------

# function to bin traits and assign colors for quant. traits
make_cols <- function(trait, palette) {
  colors <- as.character(cut(trait, breaks = 5, include.lowest = TRUE, labels=brewer.pal(5, palette)))
  return(colors)
}

# function to make vector of 1s and 0s for missing data
make_na <- function(trait) {
  colors.na <- rep(0, length(trait))
  colors.na[which(is.na(trait))] <- 1
  return(colors.na)
}

# load data ---------------------------------------------------------------

### traits
# load untransformed traits
traits <- mooreaferns::fern_traits

# add rownames to traits, drop species
rownames(traits) <- traits$species
traits$species <- NULL

# drop species with NAs for more than half of traits
trait_cutoff <- ncol(traits)*0.5
traits <- traits[rowSums(is.na(traits)) < trait_cutoff, ]

### tree
phy <- mooreaferns::fern_tree

# trim to only species with trait data
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% rownames(traits))])

# get traits in same order as tips
traits <- traits[phy$tip.label,]

# convert morphotype to binary trait
# make binary morph category: 0 is noncordate, 1 is cordate
morph_binary <- rep(0, nrow(traits))
morph_binary[which(traits$morphotype == "cordate")] <- 1
traits$morphotype <- morph_binary

# log-transform traits because they vary over a large range
# (except for dissection)
traits$stipe <- log(traits$stipe + 1)
traits$length <- log(traits$length + 1)
traits$width <- log(traits$width + 1)
traits$rhizome <- log(traits$rhizome + 1)
traits$sla <- log(traits$sla + 1)
traits$pinna <- log(traits$pinna + 1)

# then rescale (log-transform got most to range of 0-6, so standardize to 0.5-5)
traits$stipe <- rescale(traits$stipe, c(0.5,5))
traits$length <- rescale(traits$length, c(0.5,5))
traits$width <- rescale(traits$width, c(0.5,5))
traits$rhizome <- rescale(traits$rhizome, c(0.5,5))
traits$dissection <- rescale(traits$dissection, c(0.5,5))
traits$sla <- rescale(traits$sla, c(0.5,5))
traits$pinna <- rescale(traits$pinna, c(0.5,5))

# set up colors for quantitative traits -----------------------------------

# make empty list to store vectors of colors (greyscale) for plotting at tips, one for each trait
colors.list <- list()
quant.traits <- c("stipe", "rhizome", "dissection", "pinna", "sla")

palettes <- rep("Greys", length(quant.traits))
for (i in 1:length(quant.traits)) {
  colors.list[[i]] <- make_cols(traits[,colnames(traits) %in% quant.traits[i]], palettes[i])
}
names(colors.list) <- quant.traits

# make list of NAs for each trait
na.list <- list()
for (i in 1:length(quant.traits)) {
  na.list[[i]] <- make_na(traits[,colnames(traits) %in% quant.traits[i]])
}
names(na.list) <- quant.traits

# set up colors for qualitative traits ------------------------------------

# make empty list to store vectors of colors that will go on the tips
colors.qual <- list()

# define my palette of colors to choose from
# list of qualitative colors from color brewer
qualcols <- brewer.pal(9, "Set1")

# not enough colors in set1, so add paired set as well
pairedcols <- brewer.pal(8, "Paired")

# sample plot to choose colors
# plot (1:length(pairedcols), 1:length(pairedcols), pch = 16, col = pairedcols)

# list of qualitative traits
qual.traits <- c("habit", "morphotype", "gemmae", "glands", "hairs")

# use shades for present / absent
palette.qual <- c(qualcols[c(3,7)], #habit: green, brown
                  pairedcols[c(1,2)], #morphotype: dark and light blue
                  pairedcols[c(4,3)], #gemmae: dark and light green
                  pairedcols[c(6,5)], #glands: red and pink
                  pairedcols[c(8,7)]) #hairs: dark and light orange

# fill in list of tip colors for each trait
colors.qual$habit <- rep(palette.qual[1], length(traits$habit))
colors.qual$habit[traits$habit == "terrestrial"] <- palette.qual[2]
colors.qual$habit[is.na(traits$habit)] <- NA

colors.qual$morphotype <- rep(palette.qual[3], length(traits$morphotype))
colors.qual$morphotype[traits$morphotype == 0] <- palette.qual[4]
colors.qual$morphotype[is.na(traits$morphotype)] <- NA

colors.qual$gemmae <- rep(palette.qual[5], length(traits$gemmae))
colors.qual$gemmae[traits$gemmae == 0] <- palette.qual[6]
colors.qual$gemmae[is.na(traits$gemmae)] <- NA

colors.qual$glands <- rep(palette.qual[7], length(traits$glands))
colors.qual$glands[traits$glands == 0] <- palette.qual[8]
colors.qual$glands[is.na(traits$glands)] <- NA

colors.qual$hairs <- rep(palette.qual[9], length(traits$hairs))
colors.qual$hairs[traits$hairs == 0] <- palette.qual[10]
colors.qual$hairs[is.na(traits$hairs)] <- NA

# make list of NAs for each trait
na.list.qual <- list()
for (i in 1:length(qual.traits)) {
  na.list.qual[[i]] <- make_na(traits[,colnames(traits) %in% qual.traits[i]])
}
names(na.list.qual) <- qual.traits

# set up plot -------------------------------------------------------------

# set transparent outline
outline <- rgb(0,0,0,alpha=0)
# space between squares
spacer <- 7
# list of trait names to plot on top
trait.names <- c(names(colors.qual)[1], names(colors.list), names(colors.qual)[2:5])
# ymax position of legend
leg.ymax <- 120
leg.spacer <- 10

# total number of traits to plot
n.traits <- length(qual.traits) + length(quant.traits)
n.quant <- length(quant.traits)

# output plot -------------------------------------------------------------

pdf(file="traits_on_tree.pdf", width=6.5, height=9)

# base plot
plot(phy, no.margin=FALSE, cex=0.35, show.tip.label = TRUE, label.offset=(n.traits + 1)*spacer)

# first trait: add green/brown squares for epis
tiplabels(pch=22,cex=0.8, col = outline, bg = colors.qual$habit, adj = 1*spacer)

# add squares with grey heat-map for quant. traits
for (i in 1:length(colors.list)) {
  tiplabels(pch=22,cex=0.8,col = outline, bg = colors.list[[i]], adj = (i+1)*spacer)
  tiplabels(pch="\\", cex=0.5*na.list[[i]],  col = "black", adj = (i+1)*spacer)
}

# add colored squares for qual. traits
for (i in 2:length(colors.qual)) {
  tiplabels(pch=22,cex=0.8,col = outline, bg = colors.qual[[i]], adj = (i+n.quant)*spacer)
  tiplabels(pch="\\", cex=0.5*na.list.qual[[i]],  col = "black", adj = (i+n.quant)*spacer)
}

# add titles for traits
# one way: in margin, vertical (can't rotate mtext)
#for (i in 1:length(trait.names)) {
  #mtext(trait.names[i], side = 3, at = 392 + i*spacer, srt=60, las=2, cex=0.5, padj=0, adj = 0, line=-1)
#}

# another way: rotated
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
par(xpd = TRUE) #Draw outside plot area
for (i in 1:length(trait.names)) {
  text(x = 392+ i*spacer, y = corners[4]-1, trait.names[i], srt = 70, cex=0.5, adj=0)
}

# add node labels
# epiphytic clades
nodelabels ("H", getMRCA(phy, c("Hymenophyllum_multifidum", "Hymenophyllum_flabellatum")), cex=.8, bg=palette.qual[1], col="black", frame="circle")
nodelabels ("T", getMRCA(phy, c("Callistopteris_apiifolia", "Crepidomanes_bipunctatum")), cex=.8, bg=palette.qual[1], col="black", frame="circle")
nodelabels ("V", getMRCA(phy, c("Antrophyum_reticulatum", "Haplopteris_elongata")), cex=.8, bg=palette.qual[1], col="black", frame="circle")
nodelabels ("A", getMRCA(phy, c("Asplenium_affine", "Asplenium_gibberosum")), cex=.8, bg=palette.qual[1], col="black", frame="circle")
nodelabels ("E", getMRCA(phy, c("Elaphoglossum_samoense", "Elaphoglossum_savaiense")), cex=.8, bg=palette.qual[1], col="black", frame="circle")
nodelabels ("P", getMRCA(phy, c("Selliguea_plantaginea", "Oreogrammitis_raiateensis")), cex=.8, bg=palette.qual[1], col="black", frame="circle")
# others
nodelabels ("EI ", getMRCA(phy, c("Leucostegia_pallida", "Oreogrammitis_raiateensis")), cex=.8, col="black", bg="grey")
nodelabels ("EII ", getMRCA(phy, c("Diplazium_ellipticum", "Asplenium_affine")), cex=.8, col="black", bg="grey")

# add legend
legend(c(12,leg.ymax), legend = c("Epiphytic", "Terrestrial"), title = "Growth Habit", pch = 22, col = outline, fill = c(palette.qual[1], palette.qual[2]), cex = 0.7, pt.cex = 1, bty = "n", title.adj=0)
legend(c(12,leg.ymax-(1*leg.spacer)), legend = c("Cordate", "Non-cordate"), title = "Morphotype", pch = 22, col = outline, fill = c(palette.qual[3], palette.qual[4]), cex = 0.7, pt.cex = 1, bty = "n", title.adj=0)
legend(c(12,leg.ymax-(2*leg.spacer)), legend = c("Present", "Absent"), title = "Gemmae", pch = 22, col = outline, fill = c(palette.qual[5], palette.qual[6]), cex = 0.7, pt.cex = 1, bty = "n", title.adj=0)
legend(c(12,leg.ymax-(3*leg.spacer)), legend = c("Present", "Absent"), title = "Glands", pch = 22, col = outline, fill = c(palette.qual[7], palette.qual[8]), cex = 0.7, pt.cex = 1, bty = "n", title.adj=0)
legend(c(12,leg.ymax-(4*leg.spacer)), legend = c("Present", "Absent"), title = "Hairs", pch = 22, col = outline, fill = c(palette.qual[9], palette.qual[10]), cex = 0.7, pt.cex = 1, bty = "n", title.adj=0)

# add geologic scale bar using axisGeo from phyloch
data(strat2012)
axisPhylo()
axisGeo(GTS = strat2012, unit="period", cex=0.8, ages=FALSE)

dev.off()
