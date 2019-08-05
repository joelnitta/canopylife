# Data processing ----

#' Download Nitta et al 2017 Ecol Mono data zip file and 
#' extract needed data files
#'
#' @param dl_path Name of file to download the zip file to.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return Three unzipped data files:
#' - all_plots.csv: Community matrix for ferns of Moorea and Tahiti
#' including sporophytes (plot names with "_S") and gametophytes
#' (plot names with "_G").
#' - treepl_Moorea_Tahiti.tre: Dated tree for pteridophytes of Moorea and Tahiti
#' - sites.csv: Site metadata including name, latitude, longitude, and elevation (m)
#'
download_and_unzip_nitta_2017 <- function (dl_path, unzip_path, ...) {
  
  # Make sure the target directory exists
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(dl_path)))
  
  # Set url
  url <- "https://datadryad.org/bitstream/handle/10255/dryad.132050/data_and_scripts.zip?sequence=1"
  
  # Download zip file
  download.file(url, dl_path)
  
  # Unzip only needed data files to data/nitta_2017/
  unzip(dl_path, "data_and_scripts/Comm_Phylo_Analysis/data/all_plots.csv", exdir = unzip_path, junkpaths = TRUE)
  unzip(dl_path, "data_and_scripts/Comm_Phylo_Analysis/data/treepl_Moorea_Tahiti.tre", exdir = unzip_path, junkpaths = TRUE)
  unzip(dl_path, "data_and_scripts/shared_data/sites.csv", exdir = unzip_path, junkpaths = TRUE)
  unzip(dl_path, "data_and_scripts/shared_data/species.csv", exdir = unzip_path, junkpaths = TRUE)
  
}

#' Combine raw SLA (specific leaf area) data from two different sources
#'
#' @param sla_punch_data_path Path to raw data containing mass from
#' from weighing leaf punches.
#' @param sla_filmy_data_path Path to raw data containing leaf area from 
#' measuring leaf fragment images with ImageJ and masses of each fragment
#' (mostly filmy ferns, because they were too small to take punches).
#' @param species_list Vector of accepted species to include
#'
#' @return Tibble, with one measurement per row 
#' of SLA (in sq m per kg). Some individuals were measured more than
#' once.
#'
#' @examples
combine_raw_sla <- function (sla_punch_data_path, sla_filmy_data_path, species_list) {
  
  ### Read in raw data
  
  # Read in masses measured from leaf punches of known diameter.
  # - Size_of_punches is diameter in mm
  # - Number_of_punches is the number of punches used for weighing.
  # If only one punch was made this was entered as blank (NA) in the raw data.
  # - Mass is mass of punch (or punches) in g
  # - Specimen indicates Nitta collection number if associated with voucher specimen,
  # or punch sample number for a field season if no voucher
  sla_raw <- read_csv(sla_punch_data_path, col_types = "iccdcii") %>%
    clean_names
  
  # Read masses measured from leaf fragments of various shapes,
  # with area measured using ImageJ.
  # (these are species that were too small to take leaf punches)
  # Mostly filmy ferns, also includes some other random non-filmy ferns.
  # Note this data includes SLA already calculated for each specimen.
  sla_filmy <- read_csv(sla_filmy_data_path)  %>%
    clean_names
  
  #### Subset data
  # Get rid of species to exclude, only keep species in species list
  # Nitta 2543 Ptisana salicina is in error, it should be Nitta 2544.
  # Already have measurements for 2544 P. salicina, so it must have been measured twice.
  # So exclude 2544 P. salicina
  
  sla_raw <- sla_raw %>%
    mutate(notes = replace_na(notes, "none")) %>%
    filter(!str_detect(notes, "exclude")) %>%
    filter(species %in% species_list)
  
  sla_filmy <- sla_filmy %>%
    mutate(notes = replace_na(notes, "none")) %>%
    filter(!str_detect(notes, "exclude")) %>%
    filter(species %in% species_list)
  
  # Add "Nitta_" to specimen collection numbers with vouchers,
  # replace spaces with underscores
  sla_raw <- sla_raw %>%
    mutate(specimen = case_when(
      str_detect(specimen, "sample") ~ specimen,
      TRUE ~ paste0("Nitta_", specimen)
    )) %>%
    mutate(specimen = str_replace(specimen, " ", "_")) %>%
    # Exclude 2544 P. salicina
    filter(specimen != "Nitta_2544")
  
  sla_filmy <- sla_filmy %>%
    mutate(specimen = case_when(
      is.na(specimen) ~ NA_character_,
      TRUE ~ paste0("Nitta_", specimen)
    )) %>%
    mutate(specimen = str_replace(specimen, " ", "_"))
  
  #### Calculate SLA
  
  # For hole punches, use area of punch
  sla_raw <-
    sla_raw %>%
    mutate(
      # If number of punches for a sample is NA,
      # that means there was only one punch taken.
      number_of_punches = case_when(
        is.na(number_of_punches) ~ 1L,
        TRUE ~ number_of_punches
      ),
      # Fix area of punches for samples:
      # those with more than 2 punches should all be 2 mm diam
      # (as it says in notes).
      size_of_punches = case_when(
        number_of_punches > 1 ~ 2L,
        TRUE ~ size_of_punches
      ),
      # Calculate SLA based on measured area (in sq m per kg)
      # (size_of_punches is punch diameter in mm)
      sla = (((size_of_punches*0.05)^2*pi*number_of_punches)/mass)*0.1
    )
  
  # For filmy ferns, calculate SLA based on measured area (in sq m per kg)
  sla_filmy <-
    sla_filmy %>%
    mutate(
      area_cm2 = as.numeric(area_cm2),
      sla = (area_cm2/mass_g) * 0.1
    )
  
  # Combine the raw data (mulitple observations per individual)
  bind_rows(
    select(sla_raw, specimen, species, sla),
    select(sla_filmy, specimen, species, sla)
  ) %>%
    mutate(
      specimen = str_replace_all(specimen, " ", "_")
    )
  
}

#' Process raw morphological measurements data
#'
#' @param raw_morph_path Data to raw data file with morphological
#' measurements
#' @param species_list Vector of accepted species to include
#'
#' @return Tibble, with one measurement per individual per row 
#'
process_raw_cont_morph <- function (raw_morph_path, species_list) {
  read_csv(raw_morph_path, 
           na = c("?", "n/a", "NA", "N/A")) %>%
    clean_names() %>%
    # Get rid of species not in accepted list
    filter(species %in% species_list) %>%
    select(-comments, -notes, -exclude) %>%
    mutate(
      # Standardize formatting of sources.
      source = case_when(
        str_detect(source, "meas") ~ "measurement",
        TRUE ~ source
      ),
      # If specimen collection number was entered as a "1",
      # this is actually missing (i.e., no voucher specimen
      # made for that measurement).
      specimen = case_when(
        specimen == "1" ~ NA_character_,
        specimen == "MISSING_COLL_NUM" ~ NA_character_,
        TRUE ~ specimen
      )
    )
}

#' Process raw morphological qualitative trait data
#'
#' @param raw_morph_path Data to raw data file with morphological
#' trait observations
#' @param species_list Vector of accepted species to include
#'
#' @return Tibble, with one trait value per species
#'
process_raw_qual_morph <- function (raw_morph_path, species_list) {
  
  # Read in raw trait data
  qual_traits <- read_csv(raw_morph_path, na = c("?", "n/a", "NA", "N/A")) %>%
    clean_names() %>%
    select(species, habit = habit_binary, dissection, morphotype, glands, hairs, gemmae, source )
  
  # For some reason Prosaptia subnuda is missing. Add this.
  qual_traits <- bind_rows(
    qual_traits,
    tibble(species = "Prosaptia_subnuda", habit = 1)
  )
  
  # Format data sources to use names instead of number codes to minimize confusion
  # original scheme:
  #1 = lab obs.
  #2 = Nayar and Kauar 1971 "Gametophytes of homosporous ferns" The Botanical Review [@Nayar1971]
  #3 = Lloyd 1980 "Reproductive biology and gametophyte morphology of New World populations of Acrostichum aureum" AFJ [@Lloyd1980]
  #4 = field obs.
  #5 = genus/family level character (taxonomy)
  #6 = Zhang et al 2008 "Gametophyte Morphology and Development of Six Chinese Species of Pteris (Pteridaceae)" AFJ [@Zhang2008]
  #7 = Bierhorst 1967 "The gametophyte of Schizaea dichotoma" AJB [@Bierhorst1967]
  #8 = Martin et al 2006 "Efficient induction of apospory and apogamy in vitro in silver fern (Pityrogramma calomelanos L.)" Plant Cell Reports [@Martin2006]
  #9 = Tigerschiöld 1989 "Dehiscence of antheridia in thelypteroid ferns" Nordic Journal of Botany [@Tigerschiold1989]
  #10 = Tigerschiöld 1990 "Gametophytes of some Ceylonese species of Thelypteridaceae" Nordic Journal of Botany [@Tigerschiold1990]
  #11 = Atkinson 1975 "The gametophyte of five Old World thelypteroid ferns" Phytomorphology [@Atkinson1975]
  #12 = Chen et al 2014 "First insights into the evolutionary history of the Davallia repens complex" Blumea [@Chen2014a]
  qual_traits <-
    qual_traits %>%
    separate(source, c("source_1", "source_2"), ",", fill = "right") %>%
    mutate_at(vars(source_1, source_2), ~str_remove_all(., " ") %>% as.numeric) %>%
    mutate_at(vars(source_1, source_2), ~str_pad(., 2, "left", "0")) %>%
    mutate_at(
      vars(source_1, source_2),
      list(~case_when(
        . == "01" ~ "L", # lab observation
        . == "02" ~ "Nayar1971",
        . == "03" ~ "Lloyd1980",
        . == "04" ~ "F", # field observation
        . == "05" ~ "T", # based on taxonomy
        . == "06" ~ "Zhang2008",
        . == "07" ~ "Bierhorst1967",
        . == "08" ~ "Martin2006",
        . == "09" ~ "Tigerschiold1989",
        . == "10" ~ "Tigerschiold1990",
        . == "11" ~ "Atkinson1975",
        . == "12" ~ "Chen2014a",
        TRUE ~ as.character(.))
      )
    )
  
  # Make binary habit a factor
  qual_traits <-
    qual_traits %>%
    mutate(habit = factor(habit)) %>%
    mutate(habit = fct_recode(habit, terrestrial = "0", epiphytic = "1"))
  
  # keep only species in species list
  qual_traits %>%
    filter(species %in% species_list)
  
}

#' Combine trait data into trait matrix
#' (one row per species)
#'
#' @param sla_raw Raw SLA measurements
#' @param morph_cont_raw Raw continuous trait measurements
#' @param morph_qual_raw Qualitative trait observations
#'
#' @return Tibble
#' 
make_trait_matrix <- function(sla_raw, morph_cont_raw, morph_qual_raw) {
  
  # Calculate means of continuous data for each species
  morph_cont_mean <-
    morph_cont_raw %>%
    group_by(species) %>%
    summarise_if(
      is.numeric,
      ~mean(., na.rm = TRUE)
    )
  
  # Calculate grand mean SLA for each species
  sla_grand_mean <-
    sla_raw %>%
    # Means by individual first
    group_by(specimen, species) %>%
    summarize(
      sla = mean(sla, na.rm = TRUE)
    ) %>%
    # Then grand means by species
    group_by(species) %>%
    summarize(
      sla = mean(sla, na.rm = TRUE)
    ) 
  
  # Combine continuous and qualitative data
  reduce(
    list(
      sla_grand_mean,
      morph_cont_mean,
      morph_qual_raw
    ),
    full_join) %>% 
    mutate(morphotype = factor(morphotype)) %>%
    # Rename traits to simple version
    rename(
      stipe = stipe_length,
      length = frond_length,
      width = lamina_width,
      rhizome = rhizome_dia,
      pinna = pinna_pairs
    )
  
}

#' Process fern trait data for supp. info.
#'
#' @param sla_raw Dataframe; raw measurements of specific leaf area, including 
#' multiple measurments per specimen.
#' @param morph_cont_raw Dataframe; raw measurements of other continuous traits (frond length
#' and width, rhizome diameter, etc).
#' @param morph_qual_raw Dataframe; observations of qualitative traits
#' 
#' @return data frame formatted for export as CSV to be included with SI.
process_trait_data_for_si <- function (sla_raw, morph_cont_raw, morph_qual_raw) {
  
  ### SLA ###
  
  # Calculate mean SLA for each individual
  sla_mean <- 
    sla_raw %>%
    group_by(species, specimen) %>%
    summarize(sla = mean(sla, na.rm = TRUE) , sd = sd(sla, na.rm = TRUE), n = n()) %>%
    ungroup()
  
  # Calculate grand mean for each species based on individual means
  sla_grand_mean <- 
    sla_mean %>%
    group_by(species) %>%
    summarize(mean = mean(sla, na.rm = TRUE) , sd = sd(sla, na.rm = TRUE), n = n()) %>%
    mutate(trait = "sla")
  
  ### Other traits ###
  # Load raw measurements, in long format.
  # Includes only one measurment per individual but multiple individuals 
  # per species.
  morph_long <- 
    morph_cont_raw %>%
    select (-source, -specimen) %>%
    gather(key = trait, value = measurement, -species) %>%
    mutate(trait = str_replace_all(trait, "\\.", "_")) %>%
    na.omit()
  
  # Calculate means
  morph_mean <- 
    morph_long %>%
    group_by(species, trait) %>%
    summarise(mean = mean(measurement) , sd = sd(measurement), n = n())
  
  # Merge with SLA grand means
  morph_mean <- bind_rows(morph_mean, sla_grand_mean)
  
  # Format digits for printing.
  morph_mean <- 
    morph_mean %>%
    # Set number of signif digits differently for numeric vs integer measurements.
    # Numeric measurements made with normal ruler in cm, so should have to 2 digits 
    # after 0. Need to use formatC() to keep trailing zeros
    # https://kmyu.wordpress.com/2011/01/11/formatting-numbers-for-printing-in-r-rounding-and-trailing-zeroes/
    mutate(
      mean = case_when(
        trait == "pinna_pairs" ~ formatC(signif(mean, digits=2), 2, format="g"),
        TRUE ~ formatC(round(mean, digits=2), 2, format="f")
      ),
      sd = case_when(
        trait == "pinna_pairs" ~ formatC(round(sd, digits=1), 1, format="f"),
        TRUE ~ formatC(round(sd, digits=2), 2, format="f")
      )
    ) %>%
    mutate_at(vars(mean, sd), str_trim) %>%
    # Add column for mean +/- sd (n) formatted for printing
    mutate(
      print = case_when(
        n > 1 ~ glue("{mean} ± {sd} ({n})"),
        n == 1 ~ glue("{mean} ({n})")
      )
    ) %>%
    # Spread back into wide format
    select(species, trait, print) %>%
    spread(trait, print)
  
  ### Merge into final data table ###
  morph_table <- full_join(morph_mean, morph_qual_raw) %>%
    mutate(gameto_source = jntools::paste3(source_1, source_2) %>%
             str_replace(" ", ", ")) %>%
    select(-source_1, -source_2)
  
  # Add in sources for measurements
  meas_sources <- morph_cont_raw %>% 
    # Rename sources according to abbreviations
    mutate(
      source = case_when(
        source == "Ferns & Fern-Allies of South Pacific Islands" ~ "NMNS2008", 
        source == "Pteridophytes of the Society Islands" ~ "Copeland1932",
        source == "Hawaii's Ferns and Fern Allies" ~ "Palmer2003",
        source == "Flora of New South Wales" ~ "FloraNSW2019",
        source == "Brownsey 1987" ~ "Brownsey1987",
        source == "measurement" ~ "M"
      )
    ) %>%
    # Summarize by species and collapse multiple sources.
    arrange(species, source) %>%
    group_by(species) %>% 
    summarize(
      sporo_source = paste(unique(source), collapse = ", ")
    )
  
  # Add measurment sources and reformat column names/order
  left_join(morph_table, meas_sources) %>%
    select(
      species,
      `growth habit` = habit,
      morphotype,
      gemmae, 
      glands, 
      hairs, 
      `gametophyte data source` = gameto_source,
      dissection, 
      `pinna pairs` = pinna_pairs, 
      `frond length (cm)` = frond_length, 
      `frond width (cm)` = lamina_width, 
      `stipe length (cm)` = stipe_length, 
      `rhizome diam. (cm)` = rhizome_dia, 
      SLA = sla, 
      `sporophyte data source` = sporo_source
    )
  
}

#' Filter raw climate data, also calculate VPD
#'
#' @param climate_raw_path Path to raw climate data measured with
#' Hobo dataloggers (relative humidity and temperature) every 15 min
#' on Moorea and Tahiti from 2012-07-18 to 2015-02-06
#'
#' @return Tibble
process_raw_climate <- function (climate_raw_path) {
  
  # Read in raw data, includes all dataloggers from Moorea and Tahiti from
  # 2012-07-18 to 2015-02-06
  climate_raw <- read_csv(climate_raw_path)
  
  # Reformat data, filter out Tahiti sites (on Mt. Aorai)
  # Split apart date into components
  climate <- climate_raw %>%
    mutate(
      day = lubridate::day(date) %>% as.numeric,
      month = lubridate::month(date) %>% as.numeric,
      year = lubridate::year(date) %>% as.numeric) %>%
    filter(str_detect(site, "Aorai", negate = TRUE))
  
  # Split out set of days from 2013 to use for missing Tohiea 400 m ter data
  tohiea_400m_2013 <- climate %>%
    filter(date >= "2013-03-12" & date < "2013-07-06" & site == "Tohiea_400m_ter")
  
  # Replace the year 2013 with 2014 only for these dates in the Tohiea 400 m ter data
  tohiea_400m_as_2014 <-
    tohiea_400m_2013 %>%
    mutate(year = year + 1) %>%
    mutate(date = paste(year, month, day, sep = "-") %>% as.Date)
  
  # Do final filtering
  climate_filtered <- climate %>%
    # Remove bad Tohiea data from 2014
    filter(!(date >= "2014-03-12" & date < "2014-07-06" & site == "Tohiea_400m_ter")) %>%
    # Replace it with data from 2013
    bind_rows(tohiea_400m_as_2014) %>%
    # Subset to 2013-07-07 to 2014-07-05
    # - 2012-07-18 = first day have dataloggers on Tohiea
    # - 2013-07-06 = first day had all dataloggers set up on Moorea
    filter(date >= "2013-07-07" & date <= "2014-07-05") %>%
    # Delete all Rotui 600m epiphytic data 
    # (failed very early, missing for almost entire survey period)
    filter(site != "Rotui_600m_epi") %>%
    # Delete all Mouaputa 800m epiphytic data (missing from 7/10/13-8/10/13, 
    # then RH failed around 12/5/13, temp also looks suspicious)
    filter(site != "Mouaputa_800m_epi") %>%
    # Delete another chunk missing a lot of data, from 2013-11-19 to 2014-03-21
    # ('2013-11-19'), as.Date('2014-03-21')
    filter(!(date >= "2013-11-19" & date <= "2014-03-21")) %>%
    # missing data for Tohiea_800m_ter: 
    # 2014-06-15, 2014-06-16, 2014-06-19, 2014-06-20 (RH = 1%)
    filter(
      date != "2013-08-08", # RH = 1% at Tohiea_600m_ter
      date != "2014-06-15",
      date != "2014-06-16",
      date != "2014-06-19",
      date != "2014-06-20")
  
  # Reformat filtered data
  climate_filtered %>%
    mutate(
      habit = case_when(
        str_detect(site, "epi") ~ "epiphytic",
        str_detect(site,  "ter") ~ "terrestrial"
      ),
      habit = fct_relevel(habit, c("epiphytic", "terrestrial")),
      site = str_remove_all(site, "_epi|_ter"),
      vpd = plantecophys::RHtoVPD(RH, temp)
    )
  
}

#' Process raw community matrix 
#'
#' @param community_matrix_path Path to community matrix of Nitta et al 2017
#' @param species_list Vector of species to include
#' @param moorea_sites Dataframe of sites on Moorea in this study
#'
#' @return Tibble (community data matrix)
#' 
process_community_matrix <- function (community_matrix_path, species_list, moorea_sites) {
  read_csv(community_matrix_path) %>%
    rename(site = X1) %>% 
    gather(species, abundance, -site) %>%
    filter(str_detect(site, "_S")) %>%
    mutate(site = str_remove_all(site, "_S")) %>%
    # Keep only sites on Moorea
    filter(site %in% moorea_sites$site) %>%
    # Keep only species in species list (ferns of Moorea)
    filter(species %in% species_list) %>%
    spread(site, abundance)
}

#' Calculate daily means in moorea climate data
#'
#' @param moorea_climate_raw Dataframe; raw climate measurements
#'
#' @return Dataframe of daily climate measurements (daily mean, SD,
#' min, max of temperature and RH)
get_daily_means <- function (climate_raw) {
  climate_raw %>%
    group_by(site, date, habit) %>%
    summarize_at(
      vars(temp, RH, vpd),
      list(
        max = ~ max(., na.rm = TRUE),
        mean = ~ mean(., na.rm = TRUE),
        min = ~ min(., na.rm = TRUE),
        sd = ~ sd(., na.rm = TRUE)
      )
    ) %>%
    ungroup() %>%
    mutate(
      is_outlier = case_when(
        site == "Rotui_800m_slope" ~ "yes",
        TRUE ~ "no"
      )
    )
}

#' Calculate grand means of daily values for temperature and RH by site
#'
#' @param daily_means Dataframe; daily mean climate
#'
#' @return Dataframe; grand climate means of daily mean, SD,
#' min, max of temperature and RH.
get_grand_means <- function(daily_means) {
  daily_means %>%
    group_by(site, habit) %>%
    summarize_if(is.numeric, ~ mean(., na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      is_outlier = case_when(
        site == "Rotui_800m_slope" ~ "yes",
        TRUE ~ "no"
      )
    )
}

# Climate ----

#' Check for correlation between climate variables
#' and retain non-correlated ones
#'
#' @param climate_data Grand means of climate variables
#' by site
#'
#' @return Grand means of climate variables subsetted
#' to only those that are correlated less than 0.9
#' 
select_climate_vars <- function (climate_data) {
  
  # Calculate total correlation for all climatic variables
  climate_corr_all <-
    climate_data %>%
    # filter(is_outlier == "no") %>%
    select_if(is.numeric) %>%
    corrr::correlate()
  
  # Mean and min temp are correlated
  temp_corr <-
    climate_corr_all %>%
    corrr::focus(contains("temp"), mirror = TRUE)
  
  # (mean and min RH) and (min and SD RH) are correlated.
  rh_corr <-
    climate_corr_all %>%
    corrr::focus(contains("RH"), mirror = TRUE)
  
  # max VPD and sd VPD are correlated.
  vpd_corr <-
    climate_corr_all %>%
    corrr::focus(contains("vpd"), mirror = TRUE)
  
  # VPD and rel. hum. are almost all correlated.
  temp_rh_corr <-
    climate_corr_all %>%
    corrr::focus(contains("vpd"))
  
  # Drop RH entirely since VPD is a more direct measure of drying
  # intensity and almost all RH / VPD measures are correlated.
  # Drop other variables correlated > 0.9.
  climate_corr_select <-
    climate_corr_all %>%
    corrr::focus(-contains("RH"), -temp_min, -vpd_max, mirror = TRUE)
  
  # Double-check that no correlations > 0.9 remain,
  # and keep these as the final variables for analysis.
  selected_climate_vars <-
    climate_corr_select %>%
    assertr::assert(function (x) x < 0.9 | is.na(x), -starts_with("rowname")) %>%
    pull(rowname)
  
  climate_data %>%
    select(site, habit, selected_climate_vars, is_outlier)
}

# PCA ----

#' Transform traits
#'
#' @param traits dataframe; trait data with one column per trait.
#' @param log_trans logical; should log-transform be applied?
#' @param scale_traits logical; should traits be scaled?
#' @param small_number Arbitrarily small number to use in place
#' of 0 before log-transform
#' @param trans_select Character vector of trait names to log
#' transform.
#' @param scale_select Character vector of trait names to rescale.
#'
#' @return dataframe
#' 
transform_traits <- function (traits, 
                              log_trans = TRUE, 
                              scale_traits = TRUE, 
                              small_number = 0.1, 
                              trans_select = c("length", "rhizome", "stipe", "width", "pinna"), 
                              scale_select = c("sla", "dissection", "stipe", "length", 
                                               "width", "rhizome", "pinna")
) {
  
  # Log-transform
  if (log_trans == TRUE) {
    traits <-
      traits %>%
      verify(trans_select %in% colnames(traits)) %>%
      assert(is.numeric, trans_select) %>%
      # Replace zeros with arbitrarily small number
      mutate_at(trans_select, ~ifelse(. == 0, small_number, .)) %>%
      mutate_at(trans_select, log)
  }
  
  # Rescale by dividing original value by
  # its standard deviation across the dataset
  if (scale_traits == TRUE) {
    traits <-
      traits %>% 
      verify(scale_select %in% colnames(traits)) %>%
      assert(is.numeric, scale_select) %>%
      mutate_at(
        scale_select, ~scale(.) %>% as.numeric
      )
  }
  
  traits
  
}

#' Run principal components analysis on traits
#'
#' @param traits Dataframe of pre-processed traits
#' @param phy Phylogeny
#' @param cont_traits Vector of names of continuous traits
#' to include in the analysis 
#'
#' @return List of dataframes including
#'   species_locs: Positions of species along PC axes
#'   traits_locs: Positions of traits along PC axes
#'   variance: Contribution to variance of each PC
#'   
#' For each dataframe, the type of PCA is included:
#' "standard" for standard PCA (not using phylogeny), 
#' or "phylo" for PCA including phylogeny.
#' 
run_trait_PCA <- function (traits, phy, 
                           cont_traits = c("sla", "stipe", "length", "width", 
                                           "rhizome", "pinna", "dissection")) {
  
  ### Data wrangling ###
  
  # Prepare trait data
  traits <- traits %>%
    # Make sure all selected continuous traits are in the trait data
    verify(all(cont_traits %in% colnames(.))) %>%
    # Log-transform and scale
    transform_traits(
      trans_select = c("rhizome"),
      scale_select = c("sla", "dissection", "stipe", "length", 
                       "width", "rhizome", "pinna"),
      log_trans = FALSE,
      scale_traits = TRUE
    ) %>%
    # Subset to continuous traits
    select(species, habit, cont_traits) %>%
    # Keep only completely sampled species
    filter(complete.cases(.)) %>%
    # Exclude A. evecta
    filter(species != "Angiopteris_evecta") %>%
    # Keep only species in phylogeny, in phylogenetic order
    match_traits_and_tree(traits = ., phy = phy, "traits") 
  
  # Trim to only species with trait data
  phy <- match_traits_and_tree(traits, phy, "tree") 
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  # Both standard and phylo PCA expect data frames with row names
  traits_df_for_pca <- select(traits, species, cont_traits) %>%
    as.data.frame
  rownames(traits_df_for_pca) <- traits_df_for_pca$species
  traits_df_for_pca$species <- NULL
  
  ### Standard PCA ###
  pca_std_results <- PCA(traits_df_for_pca, graph=FALSE)
  
  # Get position of species along PC axes
  pca_std_species_locs <- pca_std_results$ind$coord %>%
    as.data.frame %>%
    rownames_to_column("species") %>%
    as_tibble %>%
    rename_if(is.numeric, ~str_replace_all(., "Dim\\.", "PC")) %>%
    mutate(
      analysis_type = "standard"
    )
  
  # Get position of traits along PC axes
  pca_std_trait_locs <- pca_std_results$var$coord %>%
    as.data.frame %>%
    rownames_to_column("trait") %>%
    as_tibble %>%
    rename_if(is.numeric, ~str_replace_all(., "Dim\\.", "PC")) %>%
    mutate(
      analysis_type = "standard"
    )
  
  # Get contribution to variance of each PC
  pca_std_variance <- t(pca_std_results$eig) %>%
    as.data.frame() %>%
    rownames_to_column("variance_type") %>%
    as_tibble %>%
    rename_if(is.numeric, ~str_replace_all(., "comp ", "PC")) %>%
    filter(variance_type  != "eigenvalue") %>%
    mutate(variance_type = case_when(
      variance_type == "percentage of variance" ~ "Proportion of Variance",
      variance_type == "cumulative percentage of variance" ~ "Cumulative Proportion")
    ) %>%
    mutate_if(is.numeric, ~magrittr::multiply_by(., 0.01)) %>%
    mutate(
      analysis_type = "standard"
    )
  
  ### Phylogenetic PCA ###
  pca_phy_results <- phyl.pca(phy, traits_df_for_pca, method="BM")
  
  # Get position of species along PC axes
  pca_phy_species_locs <- pca_phy_results$S %>%
    as.data.frame %>%
    rownames_to_column("species") %>%
    as_tibble %>%
    mutate(
      analysis_type = "phylogenetic"
    )
  
  # Get position of traits along PC axes
  pca_phy_trait_locs <- pca_phy_results$L %>%
    as.data.frame %>%
    rownames_to_column("trait") %>%
    as_tibble %>%
    mutate(
      analysis_type = "phylogenetic"
    )
  
  # Get contribution to variance of each PC
  pca_phy_variance <- 
    pca_phy_results %>%
    summary %>%
    magrittr::extract2("importance") %>%
    as.data.frame %>%
    rownames_to_column("variance_type") %>%
    as_tibble %>%
    mutate(
      analysis_type = "phylogenetic"
    )
  
  ### Combine results as list
  list(
    species_locs = bind_rows(pca_phy_species_locs, pca_std_species_locs),
    traits_locs = bind_rows(pca_phy_trait_locs, pca_std_trait_locs),
    variance =  bind_rows(pca_phy_variance, pca_std_variance)
  )
}

# Phylogentic signal ----

#' Analyze phylogenetic signal in a continuous trait of interest
#'
#' @param selected_trait Name of trait to analyze phylogenetic signal
#' @param traits Dataframe including all untransformed traits, with
#' 'species' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return List of estimated Blomberg's K and Pagel's lambda and
#' their significance
#' 
analyze_cont_phylosig <- function (selected_trait, traits, phy) {
  
  # Trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_select <- traits %>% select(species, selected_trait) %>%
    remove_missing(na.rm = TRUE)
  
  traits_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "traits") 
  phy_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "tree") 
  
  # Extract named vector of trait values for phylosig()
  trait_vec <- pull(traits_trim, selected_trait) %>%
    set_names(traits_trim$species)
  
  # Run phylosig() on selected trait
  # using Blomberg's K
  k_model <- phylosig(phy_trim, trait_vec, method = "K", test = TRUE)
  # and Pagel's lambda
  lambda_model <- phylosig(phy_trim, trait_vec, method = "lambda", test = TRUE)
  
  # get model results
  list(trait = selected_trait,
       kval = k_model$K,
       k_pval = k_model$P,
       lambda = lambda_model$lambda,
       lambda_pval = lambda_model$P)
  
}

#' Analyze phylogenetic signal in binary traits
#'
#' @param traits Dataframe including all untransformed traits, with
#' 'species' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return Dataframe of estimated Fritz and Purvis' D values and
#' associated statistics for each binary trait in `traits`
#' 
analyze_binary_phylosig <- function (traits, phy) {
  
  # Format binary traits as integer
  traits_binary <- 
    traits %>% 
    mutate (
      morphotype = case_when (
        morphotype == "cordate" ~ 1,
        morphotype != "cordate" ~ 0 ),
      habit = case_when (
        habit == "terrestrial" ~ 0,
        habit == "epiphytic" ~ 1)) %>%
    mutate_at(vars(habit, morphotype, gemmae, glands, hairs), as.integer) %>%
    select(species, habit, morphotype, gemmae, glands, hairs)
  
  # Nest by trait and loop phylo.d()
  # Note: some of the traits include NAs.
  # For phylo.d(), need to have matching tree and trait data with all NAs removed
  # Don't do comparative.data() on all the traits together, or species missing 
  # data for ANY trait will get dropped
  traits_binary %>%
    gather(trait, value, -species) %>%
    nest(-trait) %>%
    mutate(
      # Construct a comparative data object for each
      # binary trait
      comp_data = map(
        data, 
        ~ comparative.data(
          phy = phy, 
          data = as.data.frame(.), 
          names.col= "species", 
          na.omit=TRUE)
      ),
      # Run phylo.d on each binary trait
      phylo_d_out = map(
        comp_data,
        ~phylo.d(data = ., binvar = value) # phylo.d() uses NSE for `binvar`
      ),
      # Extract the useful bits of information from each model fit
      phylo_d_summary = map(
        phylo_d_out,
        ~tibble(
          num_present = pluck(., "StatesTable", 1), 
          num_absent = pluck(., "StatesTable", 2), 
          D = pluck(., "DEstimate"),
          prob_random = pluck(., "Pval1"), 
          prob_brownian = pluck(., "Pval0")
        )
      )
    ) %>%
    select(trait, phylo_d_summary) %>%
    unnest()
  
}

# Correlated evolution ----

#' Test for correlated evolution between growth habit and another
#' binary trait of interest
#'
#' @param selected_trait Name of trait to test for correlated evolution
#' with growth habit
#' @param traits Dataframe including all untransformed traits, with
#' 'species' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return List of loglikelihood values for the independent
#' and dependent model and p value for the test if they
#' are different.
#' 
test_corr_evo <- function (selected_trait, traits, phy) {
  
  # Make morphotype binary
  traits <- mutate (traits, morphotype = case_when (
    morphotype == "cordate" ~ "cordate",
    morphotype != "cordate" ~ "noncordate"))
  
  # Trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_select <- traits %>% select(species, habit, selected_trait) %>%
    remove_missing(na.rm = TRUE)
  
  traits_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "traits") 
  phy_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "tree") 
  
  # Format input data for phytools
  # - Extract named vector of trait values
  trait_vec <- pull(traits_trim, selected_trait) %>%
    set_names(traits_trim$species)
  # - Extract named vector of growth habit
  habit_vec <- pull(traits_trim, habit) %>%
    set_names(traits_trim$species)
  
  # run fitPagel() on selected trait vs. growth habit
  model <- fitPagel(tree = phy_trim, x = trait_vec, y = habit_vec)
  
  # get model results
  list(trait = selected_trait,
       logL_indep = model$independent.logL, 
       logL_dep = model$dependent.logL, 
       likelihood_ratio = model$lik.ratio, 
       pval = model$P)
  
}

#' Calculate phylogenetically independent contrasts between trait
#' values comparing between species that differ in (another) binary
#' trait
#'
#' @param traits Dataframe of traits, including species, a binary
#' trait (here, growth habit) and several other continuous traits 
#' @param phy Phylogeny of the species
#'
#' @return Dataframe
#' 
run_pic <- function (traits, phy) {
  
  ### Setup
  # Make sure growth habit is factor with levels set as we expect
  assert_that(isTRUE(all.equal(levels(traits$habit), c("terrestrial", "epiphytic"))))
  
  # Keep only species, habit, and quantitative traits to test
  traits <- traits %>%
    select(species, habit, stipe, length, width, rhizome, dissection, pinna, sla) %>%
    # Keep only species in phylogeny
    match_traits_and_tree(traits = ., phy = phy, "traits") 
  
  # Trim to only species with trait data
  phy <- match_traits_and_tree(traits, phy, "tree")
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  ### Run PICs using brunch
  traits %>%
    gather(trait, value, -species, -habit) %>%
    nest(-trait) %>%
    mutate(
      # Construct a comparative data object for each trait
      # so none are missing any data.
      comp_data = map(
        data, 
        ~ comparative.data(
          phy = phy, 
          data = as.data.frame(.), 
          names.col= "species", 
          na.omit=TRUE)
      ),
      # Run brunch on each trait
      model = map(
        comp_data,
        ~ brunch(value ~ habit, data = .) # phylo.d() uses NSE
      ),
      # Extract the useful bits of information from each model fit
      pic_table = map(model, caic.table),
      summary = map(model, summary),
      tidy_summary = map(summary, tidy),
      contrasts_summary = map(pic_table, 
                              ~ tibble(
                                n_contrasts = nrow(.),
                                n_pos_con = sum(.[,1] > 0)
                              )
      )
    ) %>%
    # Select only results and format as output dataframe
    select(trait, tidy_summary, contrasts_summary) %>%
    unnest
}

# Community and functional diversity ----

#' Convert community matrix including epiphytic and terrestrial species to
#' a community matrix with separate epiphytic and terretrial communities
#'
#' @param comm Community data matrix
#' @param traits Trait data including growth habit (column "habit", with values
#' as "terrestrial" or "epiphytic")
#' @param drop_zero_abun Logical; should zero-abundance species be dropped from
#' the final community data matrix?
#'
make_epi_ter_comm <- function (comm, traits, drop_zero_abun = FALSE) {
  
  ### Format community matrix ###
  # Split communities into epiphytic / terrestrial,
  # but keep together in a single community dataframe
  
  # - Add growth habit to community data
  comm <- left_join(
    comm,
    select(traits, species, habit)
  ) %>%
    # Make sure habit isn't missing for any species
    assert(not_na, habit)
  
  # - Rename columns with _E or _T at end for epiphytic or terrestrial
  # Set species abundances to 0, e.g. for terrestrial species 
  # in an epiphytic community.
  comm_by_habit <- left_join(
    comm %>%
      mutate_at(
        vars(-species, -habit), 
        ~ case_when(habit == "terrestrial" ~ 0, TRUE ~ .)) %>%
      rename_at(vars(-species, -habit), ~ paste0(., "_E")) %>%
      select(-habit),
    comm %>%
      mutate_at(
        vars(-species, -habit), 
        ~ case_when(habit == "epiphytic" ~ 0, TRUE ~ .)) %>%
      rename_at(vars(-species, -habit), ~ paste0(., "_T")) %>%
      select(-habit)
  ) %>%
    # For FD::dbFD(), input community data needs to be data.frame with
    # rows as sites (and rownames) and columns as species
    gather(site, abundance, -species) %>%
    group_by(species)
  
  if(isTRUE(drop_zero_abun)) comm_by_habit <- filter(comm_by_habit, sum(abundance) > 0)
  
  spread(comm_by_habit, species, abundance)
}

#' Analyze phylogenetic community structure in epiphytic and 
#' terrestrial communities separately
#'
#' @param comm Community matrix with sites as columns and species as rows
#' @param phy Phylogeny
#' @param traits Traits of each species, including growth habit
#' @param null_model Name of null model to use for picante ses.mpd and ses.mntd
#' functions
#' @param iterations Number of iterations to use for picante ses.mpd and ses.mntd
#' functions
#'
#' @return Dataframe; results of picante::ses.mpd and picante::ses.mntd
#' merged together.
#' 
analyze_phy_struc_by_habit <- function (comm, phy, traits, null_model = "phylogeny.pool", iterations = NULL) {
  
  ### Prepare data ###
  
  # Make sure all species are in the phylogeny
  # (but don't subset phylogeny to only species in Moorea
  # communities - we want to use the species pool including
  # both ferns of Moorea and Tahiti fro mpd, mntd).
  assert_that(all(comm$species %in% phy$tip.label))
  
  # Split communities into epiphytic / terrestrial taxa, 
  # but keep together in a single dataframe so we can
  # use whole fern community as null model.
  # 
  # Rename columns with _E or _T at end for epiphytic or terrestrial
  # Set species abundances to 0, e.g. for terrestrial species 
  # in an epiphytic community.
  # 
  # For ses.mpd() and ses.mntd(), needs to be a dataframe
  # with sites as rows and species as columns, and rownames
  # equal to site.
  comm_by_habit <- make_epi_ter_comm(comm, traits, drop_zero_abun = FALSE) %>%
    column_to_rownames("site")
  
  ### Run community structure analysis ###
  
  # Look at phylo structure of each plot individually, compared to null 
  # (standard effect size). 
  # Null model shuffles tips of phylogeny in plots, keeps freq and richness the same
  # Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) 
  # indicate phylogenetic evenness.
  # Negative SES values and low quantiles (mpd.obs.p < 0.05) indicate 
  # phylogenetic clustering.
  # SES values of 0 are for trees with species spread randomly across the tree.
  # mpd.obs.z is the standardized effect size of mpd vs. null communities 
  # (equivalent to -NRI)
  
  # Run ses.mpd, convert output to tibble
  mpd_out <- ses.mpd(
    comm_by_habit, 
    cophenetic(phy), 
    null.model = null_model, 
    abundance.weighted = TRUE, 
    runs = 999,
    iterations = iterations) %>%
    rownames_to_column("site") %>%
    as_tibble
  
  # Run ses.mntd, convert output to tibble
  mntd_out <- ses.mntd(
    comm_by_habit, 
    cophenetic(phy), 
    null.model = null_model, 
    abundance.weighted = TRUE, 
    runs = 999,
    iterations = iterations) %>%
    rownames_to_column("site") %>%
    as_tibble
  
  ### Merge results ###
  left_join(
    select(mpd_out, site, ntaxa, starts_with("mpd")),
    select(mntd_out, site, starts_with("mntd"))
    # Convert back to site and growth habit as separate columns
  ) %>%
    mutate(habit = case_when(
      str_detect(site, "_E") ~ "epiphytic",
      str_detect(site, "_T") ~ "terrestrial",
    )) %>%
    mutate(
      site = str_remove_all(site, "_E|_T")
    ) %>% 
    select(site, habit, everything())
}

#' Analyze phylogenetic community structure in epiphytic and 
#' terrestrial communities separately
#'
#' @param comm Community matrix with sites as columns and species as rows
#' @param phy Phylogeny
#' @param traits Traits of each species, including growth habit
#' @param null_model Name of null model to use for picante ses.mpd and ses.mntd
#' functions
#' @param iterations Number of iterations to use for picante ses.mpd and ses.mntd
#' functions
#' @param abundance_weighted Logical; should abundance weighting be used?
#'
#' @return Dataframe; results of picante::ses.mpd and picante::ses.mntd
#' merged together.
#' 
analyze_func_struc_by_habit <- function (
  comm, traits, 
  null_model = "phylogeny.pool", iterations = NULL, abundance_weighted = TRUE,
  traits_select =c("stipe", "length", "width", "rhizome",
                   "sla", "pinna", "dissection",
                   "morphotype", "glands", "hairs", "gemmae")) {
  
  ### Prepare data ###
  
  # Subset to only those species with trait data
  comm <-
    comm %>%
    filter(species %in% traits$species)
  
  # Split communities into epiphytic / terrestrial taxa, 
  # but keep together in a single dataframe so we can
  # use whole fern community as null model.
  # 
  # For ses.mpd() and ses.mntd(), needs to be a dataframe
  # with sites as rows and species as columns, and rownames
  # equal to site.
  comm_by_habit <- make_epi_ter_comm(comm, traits, drop_zero_abun = FALSE) %>%
    column_to_rownames("site")
  
  # Calculate Gower distance matrix
  # Input needs to be a dataframe
  # with species as rows and traits as columns, and rownames
  # equal to species.
  traits_df <- select(traits, species, traits_select) %>%
    column_to_rownames("species")
  
  # Match species order in community data
  traits_df <- traits_df[colnames(comm_by_habit), ]
  
  dist_mat <- FD::gowdis(traits_df)
  
  ### Run community structure analysis ###
  
  # Look at phylo structure of each plot individually, compared to null 
  # (standard effect size). 
  # Null model shuffles tips of phylogeny in plots, keeps freq and richness the same
  # Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) 
  # indicate phylogenetic evenness.
  # Negative SES values and low quantiles (mpd.obs.p < 0.05) indicate 
  # phylogenetic clustering.
  # SES values of 0 are for trees with species spread randomly across the tree.
  # mpd.obs.z is the standardized effect size of mpd vs. null communities 
  # (equivalent to -NRI)
  
  # Run ses.mpd, convert output to tibble
  mpd_out <- ses.mpd(
    comm_by_habit, 
    dist_mat, 
    null.model = null_model, 
    abundance.weighted = abundance_weighted, 
    runs = 999,
    iterations = iterations) %>%
    rownames_to_column("site") %>%
    as_tibble
  
  # Run ses.mntd, convert output to tibble
  mntd_out <- ses.mntd(
    comm_by_habit, 
    dist_mat, 
    null.model = null_model, 
    abundance.weighted = abundance_weighted, 
    runs = 999,
    iterations = iterations) %>%
    rownames_to_column("site") %>%
    as_tibble
  
  ### Tidy results ###
  left_join(
    select(mpd_out, site, ntaxa, starts_with("mpd")),
    select(mntd_out, site, starts_with("mntd"))
    # Convert back to site and growth habit as separate columns
  ) %>%
    mutate(habit = case_when(
      str_detect(site, "_E") ~ "epiphytic",
      str_detect(site, "_T") ~ "terrestrial",
    )) %>%
    mutate(
      site = str_remove_all(site, "_E|_T")
    ) %>% 
    select(site, habit, everything())
}

#' Calculate community-weighted means and standard deviations
#'
#' @param traits Traits of each species, including growth habit
#' @param comm Community matrix
#' @param site_data Dataframe of sites on Moorea in this study
#' @param traits_select Vector of trait names to include
#'
#' @return Tibble in long format
#' 
calculate_cwm <- function (traits, comm, site_data, 
                           traits_select = c(
                             "stipe", "length", "width", "rhizome",
                             "sla", "pinna", "dissection"
                           )) {
  
  comm %>%
    filter(species %in% traits$species) %>%
    gather(site, abundance, -species) %>%
    filter(abundance > 0) %>%
    left_join(traits) %>% 
    uncount(abundance) %>%
    left_join(site_data) %>%
    select(site, species, habit, el, traits_select) %>%
    gather(trait, value, -site, -species, -el, -habit) %>%
    group_by(habit, trait, site, el) %>%
    summarize(
      cwm = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE)
    ) %>%
    ungroup
  
}

#' Helper function to extract output from FD::dbFD()
#' 
#' FD::dbFD() outputs results in a list of named
#' vectors and dataframes. This function extracts
#' one of the named vectors as a dataframe.
#'
#' @param dbFD_output Output of FD::dbFD()
#' @param fd_metric Name of metric to extract
#'
#' @return Dataframe
#'
#' @examples
#' library(FD)
#' ex1 <- dbFD(dummy$trait, dummy$abun)
#' tidy_fd_output(ex1, "FDiv")
tidy_fd_output <- function (dbFD_output, fd_metric) {
  dbFD_output %>% 
    magrittr::extract2(fd_metric) %>%
    as.data.frame %>% 
    rownames_to_column("site") %>% 
    set_names("site", fd_metric)
}

#' Analyze functional diversity in epiphytic and 
#' terrestrial communities
#'
#' @param comm Community matrix with sites as columns and species as rows
#' @param traits Traits of each species, including growth habit
#'
#' @return Dataframe; results of FD::dbFD() in a single dataframe, including
#' site name and growth habit.
#' 
analyze_fd_by_habit <- function (traits, comm) {
  
  ### Format community matrix ###
  # Split communities into epiphytic / terrestrial,
  # but keep together in a single community dataframe so we can
  # analyze FD taking into account PC space of all species.
  comm_by_habit <- make_epi_ter_comm(comm, traits, drop_zero_abun = TRUE) %>%
    column_to_rownames("site") 
  
  ### Format traits ###
  # Subset traits to those in community, but don't transform
  # (for calculating CWMs)
  traits_sub <- 
    traits %>%
    # Keep only species in community data
    filter(species %in% colnames(comm_by_habit)) %>%
    # Keep only species name and numeric characters
    select(species, colnames(.)[map_lgl(., is.numeric)]) %>%
    # For FD::dbFD(), input trait data needs to be data.frame with
    # rows as species with rownames.
    column_to_rownames("species")
  
  ### Check order of community data and metadata ###
  assert_that(
    isTRUE(all.equal(sort(colnames(comm_by_habit)), sort(rownames(traits_sub)))),
    msg = "community and metadata don't match")
  
  traits_sub <- traits_sub[colnames(comm_by_habit), ]
  
  ### Run FD analysis ###
  
  # Set names of FD metrics to exctract
  fd_metrics <- c("FRic", "FEve", "FDiv", "FDis", "RaoQ") %>%
    set_names(.)
  
  # Calculate functional diversity on standardized trait data 
  # (set stand.x to TRUE)
  func_div <- dbFD(x = traits_sub, a = comm_by_habit, stand.x = TRUE,
                   corr = "lingoes", m = "max", 
                   w.abun = TRUE, calc.CWM = FALSE, stand.FRic = FALSE)
  
  
  ### Tidy the results ###
  map(fd_metrics, ~ tidy_fd_output(func_div, .)) %>%
    reduce(left_join, by = "site") %>%
    as_tibble %>% 
    # Convert back to site and growth habit as separate columns
    mutate(habit = case_when(
      str_detect(site, "_E") ~ "epiphytic",
      str_detect(site, "_T") ~ "terrestrial",
    )) %>%
    mutate(
      site = str_remove_all(site, "_E|_T")
    ) %>% 
    select(site, habit, everything())
  
}

#' Analyze functional diversity in epiphytic and 
#' terrestrial communities
#'
#' @param comm Community matrix with sites as columns and species as rows
#' @param traits Traits of each species, including growth habit
#'
#' @return Dataframe; results of FD::dbFD() in a single dataframe, including
#' site name and growth habit.
#' 
analyze_cwm_by_habit <- function (traits, comm) {
  
  ### Format community matrix ###
  # Split communities into epiphytic / terrestrial,
  # but keep together in a single community dataframe so we can
  # analyze FD taking into account PC space of all species.
  comm_by_habit <- make_epi_ter_comm(comm, traits, drop_zero_abun = TRUE) %>%
    column_to_rownames("site") 
  
  ### Format traits ###
  # Subset traits to those in community, but don't transform
  # (for calculating CWMs)
  traits_sub <- 
    traits %>%
    # Keep only species in community data
    filter(species %in% colnames(comm_by_habit)) %>%
    # Keep only species name and numeric characters
    select(species, colnames(.)[map_lgl(., is.numeric)]) %>%
    # For FD::dbFD(), input trait data needs to be data.frame with
    # rows as species with rownames.
    column_to_rownames("species")
  
  ### Check order of community data and metadata ###
  assert_that(
    isTRUE(all.equal(sort(colnames(comm_by_habit)), sort(rownames(traits_sub)))),
    msg = "community and metadata don't match")
  
  traits_sub <- traits_sub[colnames(comm_by_habit), ]
  
  ### Run FD analysis ###
  
  # Calculate community-weighted means on unstandardized trait data
  cwm <- dbFD(x = traits_sub, a = comm_by_habit, stand.x = FALSE,
              corr = "lingoes", m = "max", 
              w.abun = TRUE, calc.CWM = TRUE, stand.FRic = FALSE)
  
  ### Tidy the results ###
  cwm$CWM %>% 
    rownames_to_column("site") %>% 
    as_tibble %>%
    # Convert back to site and growth habit as separate columns
    mutate(habit = case_when(
      str_detect(site, "_E") ~ "epiphytic",
      str_detect(site, "_T") ~ "terrestrial",
    )) %>%
    mutate(
      site = str_remove_all(site, "_E|_T")
    ) %>% 
    select(site, habit, everything())
}

#' Check for correlation between diversity metrics
#' and retain non-correlated ones
#'
#' @param div_data All diversity metrics
#'
#' @return Diversity metrics subsetted
#' to only those that are correlated less than 0.9
#' 
select_div_metrics <- function (div_data) {
  
  # Calculate total correlation for all variables
  corr_all <-
    div_data %>%
    select(-site, -habit) %>%
    corrr::correlate()
  
  # Check for correlations > 0.9
  # (comment-out when not inspecting or will cause error)
  # corr_all %>%
  # assertr::assert(function (x) x < 0.9 | is.na(x), -starts_with("rowname"))
  
  # Only trait values are correlated:
  # stipe, length, width, and rhizome are all correlated with each other
  cwm_corr <-
    corr_all %>%
    corrr::focus(stipe, length, width, rhizome, mirror = TRUE)
  
  # Drop length, width, and rhizome.
  # Keep only stipe length since it has the best hypothesis
  # for functional significance
  div_corr_select <-
    corr_all %>%
    corrr::focus(-length, -width, -rhizome, mirror = TRUE)
  
  # Double-check that no correlations > 0.9 remain,
  # and keep these as the final variables for analysis.
  selected_div_vars <-
    div_corr_select %>%
    assertr::assert(function (x) x < 0.9 | is.na(x), -starts_with("rowname")) %>%
    pull(rowname)
  
  div_data %>%
    select(site, habit, selected_div_vars)
}


# Models ----


#' Make climate models and choose the best one using AIC
#'
#' Combines response_data and site_data datsets, runs a set of linear
#' models testing for the effect of elevation, growth habit, or their
#' interaction on each response variable. Chooses the best model
#' by lowest corrected AIC. Also check's model residuals for
#' spatial autocorrelation with Moran's I.
#'
#' @param data Dataframe containing response variables,
#' also must contain "site", "habit", "el", "long",
#' and "lat" columns.
#' @param resp_vars Vector of response variables in data
#' to use. Defaults to column names of data
#'
#' @return Dataframe.
choose_habit_elevation_models <- function (data, resp_vars) {
  
  # Reformat data into nested set by response variable
  data_long <- 
    data %>%
    select(site, habit, el, long, lat, resp_vars) %>%
    gather(var, value, -site, -habit, -el, -long, -lat)
  
  # Convert to spatial dataframe
  coordinates(data_long) <- c("long", "lat")
  proj4string(data_long) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # Can't nest a spatial dataframe, so use base split,
  # then make into tibble by hand
  data_list <- split(data_long, factor(data_long$var))
  
  # Make it into a tibble for looping
  data_tibble <-
    tibble(
      var = names(data_list),
      data = data_list
    )
  
  # Make all combinations of models including each response variable
  # by elevation only, by growth habit only, and by the interaction of
  # growth habit and elevation.
  all_models <-
    data_tibble %>%
    # Construct a linear model for each variable
    mutate(
      interaction = map(data, ~lm(value ~ el * habit, .)),
      habit_only = map(data, ~lm(value ~ habit, .)),
      el_only = map(data, ~lm(value ~ el, .)),
    ) %>%
    gather(model_type, model, -var, -data)
  
  # Calculate the AICc for each model,
  # and select model with lowest AICc
  all_models <-
    all_models %>% 
    mutate(AICc = map_dbl(model, sme::AICc)) %>%
    group_by(var) %>%
    arrange(var, AICc) %>%
    slice(1) %>%
    ungroup
  
  # For each model, make a distance matrix,
  # extract the residuals, and run Moran's I
  all_models %>%
    mutate(
      dist_mat = map(
        data, 
        ~make_dist_mat(.) %>% spdep::mat2listw(.)),
      residuals = map(model, residuals),
      moran_results = map2(
        .x = residuals,
        .y = dist_mat,
        ~spdep::moran.mc(.x, .y, 10000) %>% broom::tidy()
      )
    ) %>%
    select(var, model_type, model, AICc, moran_results)
  
}

# Extract model fits
extract_model_fits <- function (supported_models) {
  supported_models %>%
    mutate(
      fits = map(model, augment)
    ) %>%
    select(var, model_type, fits) %>%
    unnest
}

# Extract model parameter summaries (slope, intersect, etc)
extract_model_parameters <- function (supported_models) {
  supported_models %>%
    mutate(
      summary = map(model, tidy)
    ) %>%
    select(var, model_type, summary) %>%
    unnest
}

# Extract model summaries (pvalue, rsquared of model)
extract_model_summaries <- function (supported_models) {
  supported_models %>%
    mutate(
      summary = map(model, glance)
    ) %>%
    select(var, model_type, AICc, summary) %>%
    unnest
}

#' Run ANCOVA on climate data
#' 
#' The ANCOVA includes growth habit (factor with two levels) with
#' elevation as covariate
#'
#' @param data Dataframe including site names, growth habit
#' (epiphytic or terrestrial), elevatoin, and climate (response) variables
#' @param resp_vars Vector of response variables to test
#'
#' @return Dataframe
#' 
run_ancova <- function(data, resp_vars) {
  # Conduct ANCOVA
  data %>%
    select(site, habit, el, resp_vars) %>%
    gather(variable, value, -site, -habit, -el) %>%
    nest(-variable) %>%
    mutate(
      cov_model = map(data, ~aov(value ~ habit * el, .)),
      summary = map(cov_model, broom::tidy)
    ) %>%
    select(variable, summary) %>%
    unnest
}

#' Run a t-test comparing values between epiphytic and terrestrial
#' species
#'
#' @param data Input data. Must include column for "habit"
#' (growth habit, either epiphtyic or terrestrial).
#' @param resp_var Response variable in data
#'
#' @return Dataframe; summary of t-test results
#'
run_t_test_by_habit <- function (data, resp_var) {
  resp_var <- enquo(resp_var)
  ter_vals <- data %>% filter(habit == "terrestrial") %>% pull(!!resp_var)
  epi_vals <- data %>% filter(habit == "epiphytic") %>% pull(!!resp_var)
  t.test(x = ter_vals, y = epi_vals, alternative = "two.sided", na.action="na.omit") %>%
    tidy %>%
    mutate(var := !!resp_var) %>%
    select(var, everything())
}

#' Test a set of models specifically designed for the canopy life project
#'
#' @param data Dataframe including diversity metrics (MPD, MNTD, richness,
#' community-weighted means, etc), independent climate variables, and columns
#' for "site" and "habit".
#' @param resp_var Response variable of interest to fit models
#' @param indep_vars Vector of independent variables.
#'
#' @return List: output of FSSgam::fit.model.set
#' 
run_full_subset_canopy_mods <- function (data, resp_var, indep_vars) {
  
  # Subset data to the response variable and indep variables
  # of interest
  resp_var <- enquo(resp_var)
  
  diversity_data_select <- select(data, response := !!resp_var, site, habit, indep_vars) %>%
    # Make sure site and habit are factors
    mutate(
      site = as.factor(site),
      habit = as.factor(habit)) %>%
    # gam needs a data.frame, not a tibble
    as.data.frame()
  
  pred_vars_cont <- 
    diversity_data_select %>%
    dplyr::select(-response, -site, -habit) %>% colnames
  
  # Define basic model with random effect.
  # Set k = 3 to avoid over-fitting.
  basic_model <- mgcv::gam(
    response ~ s(site,bs="re"), 
    family = gaussian(), 
    data = diversity_data_select)
  
  # Make full model set
  model.set <- FSSgam::generate.model.set(
    use.dat = as.data.frame(diversity_data_select),
    test.fit = basic_model,
    pred.vars.cont = pred_vars_cont,
    pred.vars.fact = "habit",
    k = 3,
    null.terms="s(site,bs='re')")
  
  # Fit model set
  FSSgam::fit.model.set(model.set, max.models = 600, parallel = TRUE)
  
}

#' Extract importance of each indep var for a set of models
#'
#' @param fssgam_results Results of running FSSgam::fit.model.set()
#' on a list of models
#'
#' @return Tibble
#' 
get_important_vars <- function(fssgam_results) {
  
  importance_location <- list("variable.importance", "aic", "variable.weights.raw")
  
  map_df(fssgam_results, 
         ~ pluck(., !!!importance_location) %>% as.list %>% as_tibble,
         .id = "resp_var")
  
}

#' Get the best models for each dependent variable for a set of models
#'
#' Goodness of model fit assessed with AICc. The best-fitting models
#' within delta AICc of 3 are kept.
#'
#' @param fssgam_results Results of running FSSgam::fit.model.set()
#' on a list of models
#'
#' @return Tibble
#' 
get_best_fss_mods <- function(fssgam_results) {
  
  map_df(fssgam_results, 
         ~ pluck(., "mod.data.out") %>%
           as_tibble %>% filter(delta.AICc <= 3) %>% 
           arrange(delta.AICc),
         .id = "resp_var") %>%
    select(resp_var, modname, AICc, r2 = r2.vals, delta_AICc = delta.AICc) %>%
    arrange(resp_var, AICc)
  
}

# Spatial autocorrelation ----

#' Make a distance matrix for testing spatial autocorrelation
#'
#' @param data Spatial data.frame
#'
#' @return Matrix
#' 
make_dist_mat <- function(data) {
  # Make inverted distance matrix (so points far away have high values)
  dist_mat <- 1/as.matrix(dist(sp::coordinates(data)))
  # Set the diagonal (same site to same site) to zero
  diag(dist_mat) <- 0
  # Some sites have the exact same GPS points; 
  # repalce Inf values for these with 0
  dist_mat[is.infinite(dist_mat)] <- 0
  return(dist_mat)
}

#' Check for spatial autocorrelation in residuals of models
#' selected using full subsets analysis with Moran's I 
#'
#' @param model_set List of models; results of running full subsets analysis
#' as a loop over response variables
#' @param resp_var Response variable to select from models
#' @param best_mod Name of best model
#' @param site_data Dataframe with names and gps points (long, lat) 
#' of sites.
#'
#' @return Nested tibble with results of running Moran's I test
#' for spatial autocorrelation in residuals
#' 
check_moran_fss <- function (model_set, resp_var, best_mod, site_data) {
  
  # resp_var, best_mod must be character or will select wrong model (if factor)!
  resp_var <- as.character(resp_var)
  best_mod <- as.character(best_mod)
  
  # Extract best-fitting model residuals from fss results, and add long-lats.
  model_fits <- model_set[[resp_var]][["success.models"]][[best_mod]] %>%
    broom::augment() %>%
    dplyr::left_join(site_data)
  
  # Convert to spatial dataframe
  sp::coordinates(model_fits) <- c("long", "lat")
  sp::proj4string(model_fits) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # Make distance matrix
  dist_mat <- 1/as.matrix(dist(sp::coordinates(model_fits)))
  diag(dist_mat) <- 0
  # Some sites have the exact same GPS points; repalce Inf values for these with 0
  dist_mat[is.infinite(dist_mat)] <- 0
  
  # Conduct Moran's test on the residuals.
  # This will tell us if there is any remaining unexplained variation
  # in the model fit that is due to spatial autocorrelation
  # spdep::moran.test(model_fits$.resid, spdep::mat2listw(dist_mat)) %>% broom::tidy()
  spdep::moran.mc(model_fits$.resid, spdep::mat2listw(dist_mat), 10000) %>% broom::tidy()
}


# Misc ----

#' Match trait data and tree
#' 
#' Order of species in traits will be rearranged to match the
#' phylogeny.
#'
#' @param traits Dataframe of traits, with 'species' column and
#' additional columns, one for each trait
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the traits, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_traits_and_tree <- function (traits, phy, return = c("traits", "tree")) {
  
  assert_that("species" %in% colnames(traits))
  
  # Keep only species in phylogeny
  traits <- traits %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, traits$species))
  
  # Get traits in same order as tips
  traits <- left_join(
    tibble(species = phy$tip.label),
    traits
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  # Return traits or tree
  assert_that(return %in% c("tree", "traits"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (traits)
  }
  
}

#' Match community data and tree
#' 
#' Order of species in comm will be rearranged to match the
#' phylogeny.
#'
#' @param comm Community data frame, with one column for sites and
#' the rest for species.
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_comm_and_tree <- function (comm, phy, return = c("comm", "tree")) {
  
  assert_that("species" %in% colnames(comm))
  
  # Keep only species in phylogeny
  comm <- comm %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, comm$species))
  
  # Get comm in same order as tips
  comm <- left_join(
    tibble(species = phy$tip.label),
    comm
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))
  
  # Return comm or tree
  assert_that(return %in% c("tree", "comm"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (comm)
  }
  
}

# Function for rounding with trailing zeros
round_t <- function (x, digits) {
  round(x, digits) %>% sprintf(glue::glue("%.{digits}f"), .)
}

# Function for formatting p-values: round to 3 digits,
# or just say "< 0.001"
format_pval <- function (x, equals_sign = FALSE) {
  case_when(
    x < 0.001 ~ "< 0.001",
    isTRUE(equals_sign) ~ paste("=", round(x, 3) %>% as.character()),
    TRUE ~ round(x, 3) %>% as.character()
  )
}

#' Bind a list of dataframes and add an ID column with the
#' name of the dataframe.
#' 
#' Handy if you don't want to manually specify the names of
#' each individual dataframe.
#'
#' @param ... Input dataframes
#' @param id_col Value to use as name of column specifying
#' where the data came from (the name of each input dataframe)
#'
#' @examples
#' data_1 <- tibble(a = c(1,2), b = c(3,4))
#' data_2 <- tibble(a = c(5,6), b = c(7,8))
#' bind_data(data_1, data_2)
bind_data <- function (..., id_col = "dataset") {
  
  # Function to convert the input into a character
  # vector, where each item in the vector is the
  # name of the input
  # copied from
  # https://masterr.org/r/the-magic-of-substitute-and-deparse/
  foo1 <- function(a, ...) {
    arg = deparse(substitute(a))
    dots = substitute(list(...))[-1]
    c(arg, sapply(dots, deparse))
  }
  
  names = foo1(...) 
  
  data = list(...)
  
  names(data) <- names
  
  bind_rows(data, .id = id_col)
  
}

# Plotting ----

#' Make a single climate scatterplot
#'
#' @param yval Name of y variable
#' @param xval Name of x variable
#' @param ylab Label for y axis
#' @param xlab Label for x axis
#' @param single_line Logical; should a single trend line be used?
#' @param r_upper Logical; should the R-squared value be printed
#' in the upper-right? (FALSE means it will be printed in the
#' lower-right).
#' @param data Dataset used to plot points. Must include columns "el" for
#' elevation, "is_outlier" for outliers, "habit" for growth habit,
#' and response variable of interest.
#' @param fits Model fits
#' @param summaries Model summaries
#' @param habit_colors Colors to use for growth habit
#' @param alpha_val Transparency value to set for outlier points
#' not included in model
#'
#' @return ggplot object
#' @example
#' make_scatterplot(
#'   yval = "vpd_mean",
#'   data = climate_select,
#'   fits = climate_model_fits, 
#'   summaries = climate_model_summaries, 
#'   habit_colors = habit_colors
#' )
#' make_scatterplot(
#'   yval = "ntaxa",
#'   data = div_metrics_select,
#'   fits = div_el_model_fits, 
#'   summaries = div_el_model_summaries, 
#'   habit_colors = habit_colors
#' )
make_climate_scatterplot <- function (yval, xval = "el", ylab = yval, xlab = "Elevation (m)",
                              single_line = TRUE,
                              r_upper = TRUE,
                              data, 
                              fits, summaries, 
                              habit_colors, alpha_val = 0.5) {
  
  yval_sym <- sym(yval)
  xval_sym <- sym(xval)
  
  # Reformat model summaries for printing
  summaries <- summaries %>%
    mutate(r.squared = round_t(r.squared, 2))
  
  # Subset model fits and summaries to response variable of interest
  fits <- filter(fits, var == yval)
  summaries <- filter(summaries, var == yval)
  
  # Extract r2 and p values
  r2 <- pull(summaries, r.squared)
  p <- pull(summaries, p.value)
  
  # Set number of asterisks for printing with r2
  asterisk <- case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Need to wrap r2 plus asterisks in quotes for passing to 
  # ggplot::annotate with parse=TRUE so it won't error on the
  # asterisks
  r2 <- paste0('"', r2, asterisk, '"')
  
  # Extract model type: does the dependent variable depend on elevation
  # only, growth habit only, the interaction of both, or none?
  model_type <- pull(summaries, model_type)
  
  # Plot lines only for models with a significant elevation effect
  if(model_type == "el_only") {
    plot <- ggplot(data, aes(x = !!xval_sym))
    if (p < 0.05) plot <- plot + geom_line(data = fits, aes(y = .fitted))
    plot <- plot + geom_point(aes(y = !!yval_sym, alpha = is_outlier, color = habit))
  } else if(model_type == "interaction") {
    plot <- ggplot(data, aes(x = !!xval_sym, color = habit))
    if (p < 0.05) plot <- plot + geom_line(data = fits, aes(y = .fitted))
    plot <- plot + geom_point(aes(y = !!yval_sym, alpha = is_outlier))
  } else {
    plot <- ggplot(data, aes(x = !!xval_sym, color = habit)) +
      geom_point(aes(y = !!yval_sym, alpha = is_outlier))
  }
  
  # Set location of where to print R-squared
  if((model_type == "el_only" | model_type == "interaction") & p < 0.05) {
    if (isTRUE(r_upper)) {
      plot <- plot +
        annotate("text", x = Inf, y = Inf, 
                 label = glue::glue("italic(R) ^ 2 == {r2}"),
                 hjust = 1.1, vjust = 1.2,
                 parse = TRUE)
    } else {
      plot <- plot +
        annotate("text", x = Inf, y = -Inf, 
                 label = glue::glue("italic(R) ^ 2 == {r2}"),
                 hjust = 1.1, vjust = -0.2,
                 parse = TRUE)
    }
  }
  
  # Add the rest of the plot details
  plot +
    scale_color_manual(
      values = habit_colors
    ) +
    scale_alpha_manual(
      values = c("yes" = alpha_val, "no" = 1.0)
    ) +
    labs(
      y = ylab,
      x = xlab
    ) +
    standard_theme()
  
}

#' Make a single scatterplot
#'
#' @param yval Name of y variable
#' @param xval Name of x variable
#' @param ylab Label for y axis
#' @param xlab Label for x axis
#' @param single_line Logical; should a single trend line be used?
#' @param r_upper Logical; should the R-squared value be printed
#' in the upper-right? (FALSE means it will be printed in the
#' lower-right).
#' @param data Dataset used to plot points. Must include columns "el" for
#' elevation, "is_outlier" for outliers, "habit" for growth habit,
#' and response variable of interest.
#' @param fits Model fits
#' @param summaries Model summaries
#' @param habit_colors Colors to use for growth habit
#' @param alpha_val Transparency value to set for outlier points
#' not included in model
#'
#' @return ggplot object
#' @example
#' make_scatterplot(
#'   yval = "vpd_mean",
#'   data = climate_select,
#'   fits = climate_model_fits, 
#'   summaries = climate_model_summaries, 
#'   habit_colors = habit_colors
#' )
#' make_scatterplot(
#'   yval = "ntaxa",
#'   data = div_metrics_select,
#'   fits = div_el_model_fits, 
#'   summaries = div_el_model_summaries, 
#'   habit_colors = habit_colors
#' )
make_scatterplot <- function (yval, xval = "el", ylab = yval, xlab = "Elevation (m)",
                              single_line = TRUE,
                              r_upper = TRUE,
                              data, 
                              fits, summaries, 
                              habit_colors) {
  
  yval_sym <- sym(yval)
  xval_sym <- sym(xval)
  
  # Reformat model summaries for printing
  summaries <- summaries %>%
    mutate(r.squared = round_t(r.squared, 2))
  
  # Subset model fits and summaries to response variable of interest
  fits <- filter(fits, var == yval)
  summaries <- filter(summaries, var == yval)
  
  # Extract r2 and p values
  r2 <- pull(summaries, r.squared)
  p <- pull(summaries, p.value)
  
  # Set number of asterisks for printing with r2
  asterisk <- case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Need to wrap r2 plus asterisks in quotes for passing to 
  # ggplot::annotate with parse=TRUE so it won't error on the
  # asterisks
  r2 <- paste0('"', r2, asterisk, '"')
  
  # Extract model type: does the dependent variable depend on elevation
  # only, growth habit only, the interaction of both, or none?
  model_type <- pull(summaries, model_type)
  
  # Plot lines only for models with a significant elevation effect
  if(model_type == "el_only") {
    plot <- ggplot(data, aes(x = !!xval_sym))
    if (p < 0.05) plot <- plot + geom_line(data = fits, aes(y = .fitted))
    plot <- plot + geom_point(aes(y = !!yval_sym, color = habit))
  } else if(model_type == "interaction") {
    plot <- ggplot(data, aes(x = !!xval_sym, color = habit))
    if (p < 0.05) plot <- plot + geom_line(data = fits, aes(y = .fitted))
    plot <- plot + geom_point(aes(y = !!yval_sym))
  } else {
    plot <- ggplot(data, aes(x = !!xval_sym, color = habit)) +
      geom_point(aes(y = !!yval_sym))
  }
  
  # Set location of where to print R-squared
  if((model_type == "el_only" | model_type == "interaction") & p < 0.05) {
    if (isTRUE(r_upper)) {
      plot <- plot +
        annotate("text", x = Inf, y = Inf, 
                 label = glue::glue("italic(R) ^ 2 == {r2}"),
                 hjust = 1.1, vjust = 1.2,
                 parse = TRUE)
    } else {
      plot <- plot +
        annotate("text", x = Inf, y = -Inf, 
                 label = glue::glue("italic(R) ^ 2 == {r2}"),
                 hjust = 1.1, vjust = -0.2,
                 parse = TRUE)
    }
  }
  
  # Add the rest of the plot details
  plot +
    scale_color_manual(
      values = habit_colors
    ) +
    labs(
      y = ylab,
      x = xlab
    ) +
    standard_theme()
  
}

#' Make a single boxplot
#'
#' @param yval Name of y variable
#' @param xval Name of x variable
#' @param ylab Label for y axis
#' @param xlab Label for x axis
#' @param data Dataset used to plot points. Must include columns "el" for
#' elevation, "is_outlier" for outliers, "habit" for growth habit,
#' and response variable of interest.
#' @param t_summaries Results of t-test for response variable
#' between growth habit types.
#' @param habit_colors Colors to use for growth habit
#' 
#' @return ggplot object
#' @example
#' make_boxplot(
#'   yval = "ntaxa",
#'   data = div_metrics_select,
#'   summaries = div_t_test_results, 
#'   habit_colors = habit_colors
#' )
make_boxplot <- function (yval, ylab = yval, xlab = "Growth habit",
                          data, summaries, 
                          habit_colors) {
  
  yval_sym <- sym(yval)
  
  # Subset summaries to response variable of interest
  summaries <- filter(summaries, var == yval)
  
  # Extract p value
  p <- pull(summaries, p.value)
  
  # Set number of asterisks for printing
  asterisk <- case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Make basic boxplot by habit
  # Add asterisk for t-test result between the bars only if they
  # differ by habit (will print nothing if not signif.)
  ggplot(data, aes(x = habit)) +
    geom_boxplot(aes(y = !!yval_sym, color = habit)) +
    annotate(
      "text", label = asterisk, size = 5, fontface = "bold",
      x = 1.5, 
      y = Inf,
      vjust = 2.0, hjust = 0.5) +
    scale_color_manual(
      values = habit_colors
    ) +
    labs(
      y = ylab,
      x = xlab
    ) +
    standard_theme()
  
}

#' Make a multipart climate plot
#' 
#' Loops response variables over make_elevation_scatterplot().
#'
#' @param data Data needed for plotting raw points
#' @param resp_vars Vector of response variables (y-axis variable).
#' A single sub-plot will be generated for each.
#' @param fits Model fits
#' @param summaries Model summaries
#' @param habit_colors Colors to use for growth habit
#' @param alpha_val Transparency value to set for outlier points
#' not included in model
#'
#' @return ggplot object
#' 
make_climate_plot <- function (data, resp_vars,
                               fits, summaries, 
                               habit_colors, alpha_val = 0.5) {
  
  # Make tibble for looping individual plot function.
  # Each column is an argument for make_elevation_scatterplot()
  plotting_data <- tibble(
    # Use factor here so we can specify the order of the subplots
    yval = factor(resp_vars, 
                  levels = c("temp_max", "temp_mean", "temp_sd", 
                             "vpd_min", "vpd_mean", "vpd_sd")),
    ylab = case_when(
      yval == "temp_max" ~ "Max. temp. (°C)",
      yval == "temp_mean" ~ "Mean temp. (°C)",
      yval == "temp_sd" ~ "SD temp. (°C)",
      yval == "vpd_min" ~ "Min. VPD (kPa)",
      yval == "vpd_mean" ~ "Mean VPD (kPa)",
      yval == "vpd_sd" ~ "SD VPD (kPa)",
    ),
    data = list(data),
    fits = list(fits), 
    summaries = list(summaries)
  ) %>%
    # Sort into preferred order by yval, then convert back to character
    arrange(yval) %>%
    mutate(yval = as.character(yval))
  
  # Make plots as a column of the tibble
  plots_tibble <-
    plotting_data %>%
    mutate(
      plot = pmap(
        list(yval = yval, 
             ylab = ylab, 
             data = data, 
             fits = fits, 
             summaries = summaries, 
             habit_colors = list(habit_colors),
             alpha_val = alpha_val),
        make_climate_scatterplot)
    )
  
  # Can't use map with ggplot `+`, so tweak as needed in a loop.
  for(i in 1:3) {
    plots_tibble$plot[[i]] <- plots_tibble$plot[[i]] + theme(axis.title.x = element_blank())
  }
  
  # Combine plots into single output
  wrap_plots(plots_tibble$plot, ncol = 3, nrow = 2) & theme(
    legend.position = "none"# ,
    # Tweak margins to remove whitespace between plots
    # plot.margin = margin(t = 0.10, r = 0, b = 0, l = 0.10, unit = "in")
  )
}


# Helper function for labeling PCA plot
get_label <- function(axis_labels, analysis_type_select, axis_select) {
  axis_labels %>% 
    filter(analysis_type == analysis_type_select, axis == axis_select) %>% 
    pull(label)
}

make_pca_plot <- function (pca_results, habit_colors, traits) {
  
  # Make a dataframe of axis labels that include
  # amount of variance explained by each PC
  axis_labels <-
    pca_results$variance %>%
    filter(variance_type == "Proportion of Variance") %>%
    select(analysis_type, PC1, PC2) %>%
    gather(axis, proportion, -analysis_type) %>%
    mutate(
      label = glue("{axis} ({round(proportion, 3) %>% percent})"),
      axis = case_when(
        axis == "PC1" ~ "x",
        axis == "PC2" ~ "y"
      ))
  
  # Standard trait PCA 
  a <- 
    pca_results$traits_locs %>%
    filter(analysis_type == "standard") %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_point(size=3, shape=16) +
    geom_text_repel(aes(label = trait), show.legend=FALSE, segment.colour = NA) +
    scale_x_continuous(
      limits = c(-1,1)
    ) +
    scale_y_continuous(
      limits = c(-1,1)
    ) +
    labs(
      x = get_label(axis_labels, "standard", "x"),
      y = get_label(axis_labels, "standard", "y"),
      subtitle = "(a)"
    )
  
  # Phylogenetic trait PCA 
  b <- 
    pca_results$traits_locs %>%
    filter(analysis_type == "phylogenetic") %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_point(size=3, shape=16) +
    geom_text_repel(aes(label = trait), show.legend=FALSE, segment.colour = NA) +
    scale_x_continuous(
      limits = c(-1,1)
    ) +
    scale_y_continuous(
      limits = c(-1,1)
    ) +
    labs(
      x = get_label(axis_labels, "phylogenetic", "x"),
      y = get_label(axis_labels, "phylogenetic", "y"),
      subtitle = "(b)"
    )
  
  c <-
    pca_results$species_locs %>%
    filter(analysis_type == "standard") %>%
    left_join(select(traits, species, habit)) %>%
    ggplot(aes(x = PC1, y = PC2, color = habit)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_point(size=3, shape=16) +
    scale_color_manual(
      values = habit_colors
    ) +
    labs(
      x = get_label(axis_labels, "standard", "x"),
      y = get_label(axis_labels, "standard", "y"),
      subtitle = "(c)"
    )
  
  d <-
    pca_results$species_locs %>%
    filter(analysis_type == "phylogenetic") %>%
    left_join(select(traits, species, habit)) %>%
    ggplot(aes(x = PC1, y = PC2, color = habit)) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "dark gray") +
    geom_point(size=3, shape=16) +
    scale_color_manual(
      values = habit_colors
    ) +
    labs(
      x = get_label(axis_labels, "phylogenetic", "x"),
      y = get_label(axis_labels, "phylogenetic", "y"),
      subtitle = "(d)"
    )
  
  a + b + c + d + plot_layout(ncol = 2, nrow = 2) &
    standard_theme() &
    theme(legend.position = "none")
  
}

#' Combine community-weighted mean plots into final figure
#'
#' @param args Dataframe of independent and dependent variables
#' used to make scatterplots
#' @param scatterplots List of scatterplots
#' @param boxplots List of boxplots
#'
#' @return ggplot object
combine_cwm_plots <- function (scatterplots, boxplots) {
  
  # Combine scatterplots and boxplots into tibble
  plots_df <- bind_rows(
    tibble(
      resp_var = names(scatterplots),
      plot = as.list(scatterplots),
      plot_type = "scatter"
    ),
    tibble(
      resp_var = names(boxplots),
      plot = as.list(boxplots),
      plot_type = "box"
    )
  )
  
  # Sort plots by response and plot type
  # so that patchwork::wrap_plots will output them in the
  # correct order.
  plots_df <- plots_df %>%
    filter(resp_var %in% 
             c("dissection", "sla", "stipe", "pinna")) %>%
    mutate(resp_var = factor(
      resp_var, 
      levels = c("dissection", "sla", "stipe", "pinna"))) %>%
    mutate(plot_type = factor(
      plot_type, 
      levels = c("scatter", "box"))) %>%
    arrange(resp_var, plot_type)
  
  # Remove un-needed plot features.
  # Can't use ggplot `+` with mutate(), etc., so do old-fashioned loop.
  for(i in 1:nrow(plots_df)) {
    if(plots_df$resp_var[[i]] != "pinna") {
      plots_df$plot[[i]] <- plots_df$plot[[i]] + 
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
    }
    if(plots_df$plot_type[[i]] == "box") {
      plots_df$plot[[i]] <- plots_df$plot[[i]] + 
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
    }
  }
  
  # Combine plots into single output
  wrap_plots(plots_df$plot, ncol = 2) & theme(
    legend.position = "none",
    # Tweak margins to remove whitespace between plots
    plot.margin = margin(t = 0.10, r = 0, b = 0, l = 0.10, unit = "in")
  )
  
}

#' Combine functional diveristy plots into final figure
#'
#' @param scatterplots List of scatterplots
#' @param boxplots List of boxplots
#'
#' @return ggplot object
combine_comm_div_plots <- function (scatterplots, boxplots) {
  
  # Combine scatterplots and boxplots into tibble
  plots_df <- bind_rows(
    tibble(
      resp_var = names(scatterplots),
      plot = as.list(scatterplots),
      plot_type = "scatter"
    ),
    tibble(
      resp_var = names(boxplots),
      plot = as.list(boxplots),
      plot_type = "box"
    )
  )
  
  # Sort plots by response and plot type
  # so that patchwork::wrap_plots will output them in the
  # correct order.
  plots_df <- plots_df %>%
    filter(resp_var %in% 
             c("ntaxa", "mpd.obs.z", "mntd.obs.z", "FDiv", "FEve", "FRic")) %>%
    mutate(resp_var = factor(
      resp_var, 
      levels = c("ntaxa", "mpd.obs.z", "mntd.obs.z", "FDiv", "FEve", "FRic"))) %>%
    mutate(plot_type = factor(
      plot_type, 
      levels = c("scatter", "box"))) %>%
    arrange(resp_var, plot_type)
  
  # Remove un-needed plot features.
  # Can't use ggplot `+` with mutate(), etc., so do old-fashioned loop.
  for(i in 1:nrow(plots_df)) {
    if(plots_df$resp_var[[i]] != "FRic") {
      plots_df$plot[[i]] <- plots_df$plot[[i]] + 
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
    }
    if(plots_df$plot_type[[i]] == "box") {
      plots_df$plot[[i]] <- plots_df$plot[[i]] + 
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
    }
  }
  
  # Combine plots into single output
  wrap_plots(plots_df$plot, ncol = 2) & theme(
    legend.position = "none",
    # Tweak margins to remove whitespace between plots
    plot.margin = margin(t = 0.10, r = 0, b = 0, l = 0.10, unit = "in")
  )
}

#' Plot sporophyte and gametophyte traits on phylogenetic tree
#'
#' @param traits Dataframe with traits of fern sporophytes and gametophytes
#' including growth habit and others
#' @param phy Phylogeny
#'
#' @return ggplot object
#' 
plot_traits_on_tree <- function (traits, phy, ppgi) {
  
  # Subset and transform traits
  traits <-
    traits %>%
    # Remove if NA for more than half the traits
    mutate(num_na = rowSums(is.na(.))) %>%
    filter(num_na < (ncol(.) - 1) * 0.5) %>% # minus one b/c of species name
    # Convert morphotype to binary trait
    # make binary morph category: 0 is noncordate, 1 is cordate
    mutate(morphotype = case_when(
      morphotype == "cordate" ~ "1",
      TRUE ~ "0"
    )) %>%
    # Only keep species in tree
    match_traits_and_tree(phy = phy, "traits") %>%
    # Transform and rescale traits
    # SLA is already between 0 and 1, so don't log-transform.
    # Instead of scaling between min and max range, rescale between 0 and 1.
    transform_traits(trans_select = c("dissection", "stipe", "length", 
                                      "width", "rhizome", "pinna"),
                     scale_traits = FALSE) %>%
    mutate_at(vars(dissection, sla, stipe, length, width, rhizome, pinna),
              ~rescale(., c(0,1)))
  
  # For phylogeny, only keep species in traits
  phy <- match_traits_and_tree(traits, phy, "tree")
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  ### Set up color palette
  # define my palette of colors to choose from
  # list of qualitative colors from color brewer
  qualcols <- brewer.pal(9, "Set1")
  
  # not enough colors in set1, so add paired set as well
  pairedcols <- brewer.pal(8, "Paired")
  
  # use shades for present / absent
  qual_palette <- c(
    qualcols[c(3,7)], #habit: green, brown
    pairedcols[c(1,2)], #morphotype: dark and light blue
    pairedcols[c(4,3)], #gemmae: dark and light green
    pairedcols[c(6,5)], #glands: red and pink
    pairedcols[c(8,7)], #hairs: dark and light orange
    "grey95") %>%
    set_names(., c(
      "epiphytic", "terrestrial",
      "noncordate", "cordate",
      "gemmae_present", "gemmae_absent",
      "glands_present", "glands_absent",
      "hairs_present", "hairs_absent",
      "(Missing)"
    ))
  
  ### Format data: qualitative traits
  
  # Get y-axis order of tips when plotting tree
  phy_order <- ggtree(phy) %>% extract2("data") %>% select(species = label, phy_order = y)
  
  # Make tibble of qualitative traits in long format with
  # trait and species as factors (species in phy order)
  qual_traits <- traits %>%
    select(species, habit, morphotype, glands, hairs, gemmae) %>%
    mutate_all(as.character) %>%
    mutate_at(
      vars(morphotype),
      ~ str_replace(., "0", "noncordate") %>%
        str_replace(., "1", "cordate")
    ) %>%
    mutate_at(
      vars(glands, hairs, gemmae),
      ~ str_replace(., "0", "absent") %>%
        str_replace(., "1", "present")
    ) %>%
    mutate(
      glands = case_when(!is.na(glands) ~ paste0("glands_", glands)),
      hairs = case_when(!is.na(hairs) ~ paste0("hairs_", hairs)),
      gemmae = case_when(!is.na(gemmae) ~ paste0("gemmae_", gemmae))
    ) %>%
    gather(trait, value, -species) %>%
    # reorder species by phylo order when plotting tree
    left_join(phy_order) %>%
    mutate(
      species = fct_reorder(species, phy_order),
      trait = factor(trait, levels = c("habit", "morphotype", "gemmae", "glands", "hairs")),
      value = factor(value, levels = c("epiphytic", "terrestrial",
                                       "noncordate", "cordate",
                                       "gemmae_present", "gemmae_absent",
                                       "glands_present", "glands_absent",
                                       "hairs_present", "hairs_absent"
      )) %>%
        fct_explicit_na()
    )
  
  # Make tibble of quantitative traits in long format
  quant_traits <-
    traits %>%
    select(species, habit, c("sla", "dissection", "stipe", "rhizome", "pinna")) %>%
    gather(trait, value, -species, -habit) %>%
    # reorder species by phylo order when plotting tree
    left_join(phy_order) %>%
    mutate(species = fct_reorder(species, phy_order))
  
  ### Format data: phylogenetic tree
  
  # First make table of tips with taxonomic info added.
  # This will be used to make data for plotting families
  # as bars and labeled nodes in the tree.
  phy_tax_data <-
    ggtree(phy) %>% extract2("data") %>%
    rename(species = label) %>%
    mutate(genus = str_split(species, "_") %>% map_chr(1)) %>%
    left_join(
      select(ppgi, genus, family)
    ) %>%
    mutate(
      family = case_when(
        genus == "Amphineuron" ~ "Thelypteridaceae",
        genus == "Wibelia" ~ "Davalliaceae",
        genus == "Humata" ~ "Davalliaceae",
        genus == "Belvisia" ~ "Polypodiaceae",
        TRUE ~ family
      )
    )
  
  # Make tibble of fern families to plot along y-axis.
  # Order must match tips of tree.
  family_bars_data <-
    phy_tax_data %>%
    remove_missing() %>%
    group_by(family) %>%
    summarize(
      start = min(y),
      end = max(y),
      count = n()
    ) %>%
    mutate(
      size_class = case_when(
        count > 5 ~ 3,
        count > 2 ~ 2,
        TRUE ~ 1
      )
    ) %>%
    ungroup() %>%
    mutate(
      family = fct_reorder(family, end)
    ) %>%
    rowwise() %>%
    mutate(
      mid = mean(c(start, end))
    ) %>%
    arrange(family) %>%
    ungroup() %>%
    mutate(
      rownum = 1:nrow(.),
      is_odd = case_when(
        (rownum %% 2) != 0 ~ "yes",
        (rownum %% 2) == 0 ~ "no"
      )) %>%
    select(-rownum) %>%
    mutate(
      start = start-0.5,
      end = end+0.5)
  
  ### Subplots 1: Heatmap plot of trait values
  
  qual_heatmap <-
    ggplot(qual_traits, aes(trait, species)) +
    geom_tile(aes(fill = value)) +
    scale_fill_manual(
      values = qual_palette,
      breaks = names(qual_palette),
      labels = names(qual_palette) %>%
        str_to_sentence() %>%
        str_replace_all("_", " ")) +
    scale_x_discrete(
      position = "top",
      breaks = qual_traits$trait,
      labels = str_to_sentence(qual_traits$trait)) +
    guides(fill = guide_legend(
      title = "Traits",
      nrow = 2,
      title.position = "top")) +
    jntools::blank_y_theme() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 24/.pt),
      axis.title.x = element_blank(),
      legend.text = element_text(size = 20/.pt),
      legend.title = element_text(size = 20/.pt)
    )
  
  # Don't use frond length or width, as these are correlated with stipe length
  quant_heatmap <-
    ggplot(quant_traits) +
    geom_tile(aes(trait, species, fill = value)) +
    scale_y_discrete(
      position = "right",
      labels = function(x) str_replace_all(x, "_", " "),
      expand = c(0,0)) +
    # jntools::blank_y_theme() +
    scale_x_discrete(
      position = "top",
      breaks= c("sla", "dissection", "stipe", "rhizome", "pinna"),
      labels= c("SLA", "Dissection", "Stipe length", "Rhizome diam.", "Pinna num.")) +
    scale_fill_viridis(
      option = "B",
      breaks = c(0, 1),
      labels = c("Low", "High")) +
    labs(
      y = "",
      fill = "Relative Trait Value") +
    guides(fill = guide_colorbar(title.position = "top")) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 24/.pt),
      axis.title.x = element_blank(),
      axis.text.y = element_text(face = "italic", size = 12/.pt),
      legend.text = element_text(size = 20/.pt),
      legend.title = element_text(size = 20/.pt)
    )
  
  # Extract legend
  quant_legend <- cowplot::get_legend(quant_heatmap)
  quant_heatmap <- quant_heatmap + theme(legend.position = "none")
  
  ### Subplot 2: Phylogenetic tree
  
  # Make named vector of nodes to label
  # node is the most common recent ancestor of
  # the group we want to label, and the name
  # of the vector is the text to print for that node.
  nodes_to_label <-
    c(
      A = phy_tax_data %>%
        filter(family == "Aspleniaceae") %>%
        pull(species) %>%
        getMRCA(phy, .),
      
      E = phy_tax_data %>%
        filter(genus == "Elaphoglossum") %>%
        pull(species) %>%
        getMRCA(phy, .),
      
      H = phy_tax_data %>%
        filter(genus == "Hymenophyllum") %>%
        pull(species) %>%
        getMRCA(phy, .),
      
      P = phy_tax_data %>%
        filter(family == "Polypodiaceae") %>%
        pull(species) %>%
        getMRCA(phy, .),
      
      V = phy_tax_data %>%
        filter(genus %in% c("Haplopteris", "Antrophyum")) %>%
        pull(species) %>%
        getMRCA(phy, .)
      
    ) %>% sort
  
  # Convert named vector into dataframe based on
  # original tip data
  nodes_to_label_df <-
    phy_tax_data %>%
    filter(node %in% nodes_to_label) %>%
    mutate(species = names(nodes_to_label))
  
  # Plot tree
  phy_plot <- ggtree(phy) +
    # Expand y scale so the tips line up with the heatmap
    scale_y_continuous(expand = c(0,0),
                       # KEY: Extend from 0.5 to 0.5 past the number of tips
                       limits = c(0.5, 0.5 + length(phy$tip.label))) +
    # Expand x scale so there's no empty space next to heatmap
    scale_x_continuous(expand = c(0,0)) +
    # Probably not necessary, but if we we're using tip labels
    # would include so they don't get chopped off past plotting region.
    coord_cartesian(clip = "off") +
    # Label epiphytic clades with green labels.
    geom_label(
      data = nodes_to_label_df,
      aes(label = species),
      fill = qual_palette[["epiphytic"]],
      size = 2,
      label.padding = unit(0.1, "lines")) +
    geom_treescale(width = 50, offset = 1, x = 0 , y = 50)
  
  ### Assemble subplots into final plot
  phy_plot + qual_heatmap + quant_heatmap + plot_layout(nrow = 1, widths = c(4,2,2)) &
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      plot.margin = margin(t = 0, r = 0, b = 0, l = -0.025, unit = "in")
    )
  
}

#' Make a heatmap of importance scores
#'
#' @param important_div_vars Importance scores output from 
#' FSSgam::fit.model.set()
#'
#' @return GGplot object
#' 
make_heatmap <- function(important_div_vars) {
  
  # Reformat data for plotting
  plot_data <-
    important_div_vars %>%
    gather(indep_var, value, -resp_var) %>%
    # Specify levels of indep and resp vars so they show up along
    # the axes in the right order
    mutate(
      indep_var = factor(
        indep_var, 
        levels = c("habit", "temp_max", "temp_mean", 
                   "temp_sd", "vpd_mean", "vpd_min", "vpd_sd")),
      resp_var = factor(
        resp_var, 
        levels = c("ntaxa", "mpd.obs.z.phy", "mntd.obs.z.phy", 
                   "mpd.obs.z.func", "mntd.obs.z.func", 
                   "dissection", "sla", "stipe", "pinna"))
    ) %>%
    # Reformat the names of the indep and resp vars so they
    # look pretty
    mutate(
      indep_var = lvls_revalue(
        indep_var, 
        new_levels = c("Habit", "Max. temp.", "Mean temp.", 
                       "SD Temp.", "Mean VPD", "Min. VPD", "SD VPD")),
      resp_var = lvls_revalue(
        resp_var, 
        new_levels = c("Richness", "MPDphy", "MNTDphy", 
                       "MPDfunc", "MNTDfunc", 
                       "Dissection", "SLA", "Stipe length", "Pinna no.")) %>%
        fct_rev()
    )
  
  # Make heatmap
  ggplot(plot_data, aes(x = indep_var, y = resp_var, fill = value)) +
    geom_tile() +
    scale_fill_scico(
      palette = "bilbao", 
      breaks = c(0, 0.25, 0.5, 0.75, 1.0), 
      limits = c(0, 1.0)) +
    labs(
      x = "Predictor",
      y = "Response",
      fill = "Importance"
    ) +
    jntools::standard_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}
