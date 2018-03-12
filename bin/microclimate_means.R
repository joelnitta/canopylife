# microclimate_means.R
# calculate mean microclimate variables for sites on Moorea

# load packages
library(mooreaferns) # moorea_climate

# set working directory
setwd(here::here())

# calculate daily maximum, mean, minimum, and SD of temperature and RH by site
daily_values <- plyr::ddply(moorea_climate, plyr::.(site, date), summarize,
                            max_temp = max(temp),
                            mean_temp = mean(temp),
                            min_temp = min(temp),
                            sd_temp = sd(temp),
                            max_RH = max(RH),
                            mean_RH = mean(RH),
                            min_RH = min(RH),
                            sd_RH = sd(RH))

# calculate grand means of daily values for temperature and RH by site
grand_means <- plyr::ddply(daily_values, "site", summarize,
                           max_temp = mean(max_temp),
                           mean_temp = mean(mean_temp),
                           min_temp = mean(min_temp),
                           sd_temp = mean(sd_temp),
                           max_RH = mean(max_RH),
                           mean_RH = mean(mean_RH),
                           min_RH = mean(min_RH),
                           sd_RH = mean(sd_RH))

### reformat data to be in site x variable format
# add a column for growth habit: daily means
daily_values$habit <- ""
daily_values$habit[grep("epi", daily_values$site)] <- "epiphytic"
daily_values$habit[grep("ter", daily_values$site)] <- "terrestrial"
daily_values$site <- gsub("_epi", "", daily_values$site)
daily_values$site <- gsub("_ter", "", daily_values$site)

# grand means
grand_means$habit <- ""
grand_means$habit[grep("epi", grand_means$site)] <- "epiphytic"
grand_means$habit[grep("ter", grand_means$site)] <- "terrestrial"
grand_means$site <- gsub("_epi", "", grand_means$site)
grand_means$site <- gsub("_ter", "", grand_means$site)
