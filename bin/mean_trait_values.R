library(tidyverse)
library(mooreaferns)

setwd(here::here())

###########
### SLA ###
###########

# load raw measurements, including multiple measurements per individual
sla.raw <- mooreaferns::sla.raw

# Nitta 2543 Ptisana salicina is in ERROR, should be Nitta 2544.
# Already have measurements for 2544 P. salicina, so it must have been measured twice.
# Exclude 2544 P. salicina
sla.raw <- sla.raw[!(sla.raw$specimen=="Nitta_2544" & sla.raw$species=="Ptisana_salicina"), ]

# calculate mean SLA for each individual
sla.mean <- 
  sla.raw %>%
  group_by(specimen) %>%
  summarize(sla = mean(sla) , sd = sd(sla), n = n()) %>%
  left_join(unique(select(sla.raw, species, specimen)))

# calculate grand mean for each species based on individual means
sla.grand.mean <- 
  sla.mean %>%
  group_by(species) %>%
  summarize(mean = mean(sla) , sd = sd(sla), n = n())
  
sla.grand.mean$trait <- "sla"

####################
### other traits ###
####################

# load raw measurements, in long format
# includes only one measurment per individual but multiple individuals per species
morph.long <- 
  mooreaferns::morph.raw %>%
  as_tibble %>%
  select (-source, -specimen) %>%
  gather(key = trait, value = measurement, -species) %>%
  na.omit()

# calculate means
morph.mean <- 
  morph.long %>%
  group_by(species, trait) %>%
  summarise(mean = mean(measurement) , sd = sd(measurement), n = n())
rm(morph.long)

# merge with sla grand means
morph.mean <- bind_rows(morph.mean, sla.grand.mean)

# split means into numeric and integer (only number of pinna pairs) measurements
morph.numeric <- filter(morph.mean, trait == "frond.length" | trait == "lamina.width" | trait == "rhizome.dia" | trait == "stipe.length" | trait == "sla") 
morph.int <- filter(morph.mean, trait == "pinna.pairs")

# set number of signif digits differently for each
# numeric measurements made with normal ruler in cm, so should have to 2 digits after 0
# need to use formatC to keep trailing zeros
# see https://kmyu.wordpress.com/2011/01/11/formatting-numbers-for-printing-in-r-rounding-and-trailing-zeroes/
morph.numeric$mean <- formatC(round(morph.numeric$mean, digits=2), 2, format="f")
morph.numeric$sd <- formatC(round(morph.numeric$sd, digits=2), 2, format="f")
# for pinna pairs (integer), keep 2 digits for mean, 1 for SD
morph.int$mean <- formatC(signif(morph.int$mean, digits=2), 2, format="g")
morph.int$sd <- formatC(round(morph.int$sd, digits=1), 1, format="f")

# combine back togeter
morph.mean <- rbind(morph.numeric, morph.int)
rm(morph.numeric, morph.int)

# add column for mean +/- sd (n) formatted for printing
morph.mean$print <- ""
morph.mean$print[morph.mean$n > 1] <- paste(morph.mean$mean[morph.mean$n > 1], " ± ", morph.mean$sd[morph.mean$n > 1], " (", morph.mean$n[morph.mean$n > 1], ")", sep="" )
morph.mean$print[morph.mean$n == 1] <- paste(morph.mean$mean[morph.mean$n == 1], " (", morph.mean$n[morph.mean$n == 1], ")", sep="")

# spread back into wide format
morph.mean <- 
  select(morph.mean, species, trait, print) %>%
  spread(trait, print)

###################################
### merge into final data table ###
###################################

morph.results <- 
  select(mooreaferns::fern_traits, species, habit, dissection, morphotype, glands, hairs, gemmae) %>%
  left_join(morph.mean)

# need to add lit sources, output as CSV




# need to decide what to do with Hymenophyllum braitwaitei (right now in SLA but not other traits)
# morph.mean <- merge(morph.mean, sla.grand.mean, by = "species", all.x=TRUE, all.y=TRUE)

# format for supp. info table
# rownames(morph.mean) <- morph.mean$species
# morph.mean$species <- NULL
