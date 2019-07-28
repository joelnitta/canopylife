# Traits ----

#' Process fern trait data for supp. info.
#'
#' @param sla_raw Dataframe; raw measurements of specific leaf area, including 
#' multiple measurments per specimen.
#' @param morph_raw Dataframe; raw measurements of other traits (frond length
#' and width, rhizome diameter, etc).
#' @param fern_traits Dataframe; pre-processed trait data including one value 
#' per species.
#' @return data frame formatted for export as CSV to be included with SI.
process_trait_data_for_si <- function (sla_raw, morph_raw, fern_traits) {
  
  ### SLA ###
  # Nitta 2543 Ptisana salicina is in error, it should be Nitta 2544.
  # Already have measurements for 2544 P. salicina, so it must have been measured twice.
  # So exclude 2544 P. salicina
  
  # Calculate mean SLA for each individual
  sla_mean <- 
    sla_raw %>%
    filter(specimen != "Nitta_2544") %>%
    group_by(species, specimen) %>%
    summarize(sla = mean(sla) , sd = sd(sla), n = n()) %>%
    ungroup()
  
  # Calculate grand mean for each species based on individual means
  sla_grand_mean <- 
    sla_mean %>%
    group_by(species) %>%
    summarize(mean = mean(sla) , sd = sd(sla), n = n()) %>%
    mutate(trait = "sla")
  
  ### Other traits ###
  # Load raw measurements, in long format.
  # Includes only one measurment per individual but multiple individuals 
  # per species.
  morph_long <- 
    morph_raw %>%
    as_tibble %>%
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
  morph_table <- 
    fern_traits %>%
    select(species, habit, dissection, morphotype, glands, hairs, gemmae, gameto_source) %>%
    left_join(morph_mean)
  
  # Add in sources for measurements
  meas_sources <- morph_raw %>% 
    # Rename sources according to abbreviations
    mutate(
      source = case_when(
        source == "Ferns & Fern-Allies of South Pacific Islands" ~ "1",
        source == "Pteridophytes of the Society Islands" ~ "2",
        source == "Hawaii's Ferns and Fern Allies" ~ "3",
        source == "Flora of New South Wales" ~ "4",
        source == "Brownsey 1987" ~ "5",
        source == "measurement" ~ "M"
      )
    ) %>%
    # Summarize by species and collapse multiple sources.
    arrange(species, source) %>%
    group_by(species) %>% 
    summarize(
      source = paste(unique(source), collapse = ", ")
    )
  
  # Add measurment sources and reformat column names/order
  left_join(morph_table, meas_sources) %>%
    select(
      species,
      `growth habit` = habit,
      dissection,
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
      `sporophyte data source` = source
    )
  
}

# Climate ----

#' Reformat Moorea climate data
#' 
#' @param moorea_climate_raw Dataframe; raw values with growth habit
#' (terrestrial or epiphytic) coded in site name
#' @return Dataframe; raw values with habit as a separate column.
reformat_raw_climate <- function (moorea_climate_raw) {
  moorea_climate_raw %>%
    # Split out habit as separate column and set levels
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
                              trans_select = c("dissection", "stipe", "length", "width", 
                                               "rhizome", "pinna"), 
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
  
  # Rescale by dividing original value by range of 
  # that value (max - min) across the dataset
  if (scale_traits == TRUE) {
    traits <-
      traits %>% 
      verify(scale_select %in% colnames(traits)) %>%
      assert(is.numeric, scale_select) %>%
      mutate_at(
        scale_select, ~ . / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
      )
  }
  
  traits
  
}

#' Run principal components analysis on traits
#'
#' @param traits Dataframe of pre-processed traits
#' @param phy Phylogeny
#' @param analysis 
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
run_trait_PCA <- function (traits, phy) {
  
  ### Data wrangling ###
  
  # Make vector of continuous traits to include in PCA
  cont_vars <- c(
    "sla", "stipe", "length", "width", "rhizome", "pinna"
  )
  
  # Prepare trait data
  traits <- traits %>%
    # Log-transform and scale
    transform_traits() %>%
    # Subset to continuous traits
    select(species, habit, sla, stipe, length, width, rhizome, pinna) %>%
    # Keep only completely sampled species
    filter(complete.cases(.)) %>%
    # Keep only species in phylogeny
    match_traits_and_tree(traits = ., phy = phy, "traits") 
  
  # Trim to only species with trait data
  phy <- match_traits_and_tree(traits, phy, "tree") 
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  # Both standard and phylo PCA expect data frames with row names
  traits_df_for_pca <- select(traits, species, cont_vars) %>%
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

#' Analyze phylogenetic community structure in epiphytic and 
#' terrestrial communities separately
#'
#' @param comm Community matrix with sites as columns and species as rows
#' @param phy Phylogeny
#' @param traits Traits of each species, including growth habit
#'
#' @return Dataframe; results of picante::ses.mpd and picante::ses.mntd
#' merged together.
#' 
analyze_comm_struc_by_habit <- function (comm, phy, traits) {
  
  ### Prepare data ###
  
  # Keep only species in community with trait data
  comm <- comm %>%
    filter(species %in% traits$species)
  
  # Add growth habit to community data
  comm <- left_join(
    comm,
    select(traits, species, habit)
    ) %>%
  assert(not_na, habit)
  
  # check for mismatches/missing species between phy and comm
  phy <- match_comm_and_tree(comm, phy, "tree")
  comm <- match_comm_and_tree(comm, phy, "comm")
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))
  
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
    gather(site, abundance, -species) %>%
    spread(species, abundance) %>%
    as.data.frame()
  
  rownames(comm_by_habit) <- comm_by_habit$site
  comm_by_habit$site <- NULL
  
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
    null.model = "phylogeny.pool", 
    abundance.weighted = TRUE, 
    runs = 999, iterations = 1000) %>%
    rownames_to_column("site") %>%
    as_tibble
  
  # Run ses.mntd, convert output to tibble
  mntd_out <- ses.mntd(
    comm_by_habit, 
    cophenetic(phy), 
    null.model = "phylogeny.pool", 
    abundance.weighted = TRUE, 
    runs = 999, iterations = 1000) %>%
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
#' terrestrial communities separately
#'
#' @param comm Community matrix with sites as columns and species as rows
#' @param traits Traits of each species, including growth habit
#'
#' @return Dataframe; results of FD::dbFD() in a single dataframe, including
#' site name and growth habit.
#' 
analyze_fd_by_habit <- function (traits, comm, habit_type = c("epiphytic", "terrestrial")) {
  
  # Only include species with trait data and 
  # with at least one occurrence in all plots.
  
  # For FD::dbFD(), input community data needs to be data.frame with
  # rows as sites (and rownames) and columns as species
  comm <-
    comm %>%
    filter(species %in% traits$species) %>%
    gather(plot, abundance, -species) %>%
    group_by(species) %>%
    filter(sum(abundance) > 0) %>%
    spread(species, abundance) %>%
    column_to_rownames("plot")
  
  # Make dataframe of traits subsetted to just
  # epiphytic or terrestrial species (specify with habit_type).
  
  # For FD::dbFD(), input trait data needs to be data.frame with
  # rows as species with rownames.
  traits_sub <- 
    traits %>%
    # Keep only species in community data
    filter(species %in% colnames(comm)) %>%
    # Subset to epiphytes
    filter(habit == habit_type) %>%
    # Keep only species name and numeric characters
    select(species, colnames(.)[map_lgl(., is.numeric)]) %>%
    # Remove if NA for more than half the traits
    mutate(num_na = rowSums(is.na(.))) %>%
    filter(num_na < (ncol(.) - 1) * 0.5) %>% # minus one b/c of species name
    column_to_rownames("species")
  
  # Make community of epiphytic or terrestrial species only
  comm_sub <- select(comm, rownames(traits_sub))
  
  # Check order of community data and metadata
  if (!(all.equal(colnames(comm_sub), rownames(traits_sub)))) {
    stop ("community and metadata don't match")
  }
  
  # Set names of FD metrics to exctract
  fd_metrics <- c("FRic", "FEve", "FDiv", "FDis", "RaoQ") %>%
    set_names(.)
  
  # Run FD analysis
  # Transform traits using default settings for func div metrics
  traits_sub_trans <- rownames_to_column(traits_sub, "species") %>%
    transform_traits %>%
    column_to_rownames("species")
  
  # Calculate functional diversity on transformed trait data
  func_div <- dbFD(x = traits_sub_trans, a = comm_sub, 
                       corr = "lingoes", m = "max", 
                       w.abun = TRUE, calc.CWM = FALSE, stand.FRic = FALSE)
  
  # Calculate community-weighted means of untransformed trait data
  cwm <- dbFD(x = traits_sub, a = comm_sub, 
                  corr = "lingoes", m = "max", 
                  w.abun = TRUE, calc.CWM = TRUE, stand.FRic = FALSE)
  
  # Tidy results, add habitat type
    map(fd_metrics, ~ tidy_fd_output(func_div, .)) %>%
    reduce(left_join, by = "site") %>%
    left_join(
      cwm$CWM %>% rownames_to_column("site") %>% select(-num_na)
    ) %>%
    mutate(habit = habit_type) %>%
    as_tibble
  
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
#' by lowest corrected AIC.
#'
#' @param data Dataframe containing response variables,
#' also must contain "site", "habit", and "el" columns.
#' @param resp_vars Vector of response variables in data
#' to use. Defaults to column names of data
#'
#' @return Dataframe.
choose_habit_elevation_models <- function (data, resp_vars = colnames(data)) {
  
  # Reformat data into nested set by response variable
  data <- 
    data %>%
    select(site, habit, el, resp_vars) %>%
    gather(var, value, -site, -habit, -el) %>%
    nest(-var)
  
  # Make all combinations of models including each response variable
  # by elevation only, by growth habit only, and by the interaction of
  # growth habit and elevation.
  all_models <-
    data %>%
    # Construct a linear model for each variable
    mutate(
      interaction = map(data, ~lm(value ~ el * habit, .)),
      habit_only = map(data, ~lm(value ~ habit, .)),
      el_only = map(data, ~lm(value ~ el, .))
    ) %>%
    select(-data) %>%
    gather(model_type, model, -var)
  
  # Calculate the AICc for each model,
  # and select model with lowest AICc
  all_models %>% 
    mutate(AICc = map_dbl(model, sme::AICc)) %>%
    group_by(var) %>%
    arrange(var, AICc) %>%
    slice(1) %>%
    ungroup
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

# Plotting ----

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
        make_scatterplot)
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
    # Don't need gameto trait data source
    select(-gameto_source) %>%
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
        levels = c("ntaxa", "mpd.obs.z", "mntd.obs.z", 
                   "FDiv", "FEve", "FRic",
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
        new_levels = c("Richness", "MPD", "MNTD", 
                       "FDiv", "FEve", "FRic",
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
