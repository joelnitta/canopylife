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
      site = str_remove_all(site, "_epi|_ter")
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
    summarize(
      max_temp = max(temp),
      mean_temp = mean(temp),
      min_temp = min(temp),
      sd_temp = sd(temp),
      max_RH = max(RH),
      mean_RH = mean(RH),
      min_RH = min(RH),
      sd_RH = sd(RH)
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
    summarize(
      max_temp = mean(max_temp),
      mean_temp = mean(mean_temp),
      min_temp = mean(min_temp),
      sd_temp = mean(sd_temp),
      max_RH = mean(max_RH),
      mean_RH = mean(mean_RH),
      min_RH = mean(min_RH),
      sd_RH = mean(sd_RH)
    ) %>%
    ungroup() %>%
    mutate(
      is_outlier = case_when(
        site == "Rotui_800m_slope" ~ "yes",
        TRUE ~ "no"
      )
    )
}

#' Make climate models
#'
#' @param climate_data 
#' @param site_data 
#'
#' @return Dataframe.
make_climate_models <- function (climate_data, site_data) {
  
  # Make dataset for modeling. Each row is a nested dataset for a single
  # climate variable so it will be easier to loop over them.
  data_for_modeling <- 
    # Combine grand climate means and site data into 'diversity data'
    left_join(climate_data, site_data) %>%
    # REMOVE OUTLIER: Mt Rotui exposed slope
    filter(is_outlier == "no") %>%
    # Reformat data into nested set by variable (mean and SD of temp and RH)
    select(site, habit, el, mean_temp, sd_temp, min_RH, sd_RH) %>%
    gather(var, value, -site, -habit, -el) %>%
    nest(-var)
  
  # Make all combinations of models including each climate variable
  # by elevation only, by growth habit only, and by the interaction of
  # growth habit and elevation.
  all_models <-
    data_for_modeling %>%
    # Construct a linear model for each variable
    mutate(
      interaction = map(data, ~lm(value ~ el * habit, .)),
      habit_only = map(data, ~lm(value ~ habit, .)),
      el_only = map(data, ~lm(value ~ el, .))
    ) %>%
    select(-data) %>%
    gather(model_type, model, -var)
  
  # Run analysis of variance on each interaction model and filter to only
  # the significant parts of the model. Use this to determine which of the
  # possible set of models is best.
  aov_results <-
    data_for_modeling %>%
    # Construct a linear model for each variable, and pull out summaries
    mutate(
      model = map(data, ~aov(value ~ el * habit, .)),
      summary = map(model, tidy)
    ) %>%
    select(var, summary) %>%
    unnest() %>%
    # Filter to only features of model that are significant
    filter(p.value < 0.05) %>%
    select(var, term) %>%
    group_by(var) %>%
    summarize(
      term = paste(unique(term), collapse = " * ")
    ) %>%
    mutate(
      model_type = case_when(
        term == "el" ~ "el_only",
        term == "el * habit" ~ "interaction",
        term == "habit" ~ "habit_only"
      )
    ) %>%
    select(var, model_type)
  
  # Based on the AOV results, filter to only the best models.
  left_join(aov_results, all_models) 
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

# Extract model summaries
extract_model_summaries <- function (supported_models) {
  supported_models %>%
    mutate(
      summary = map(model, tidy)
    ) %>%
    select(var, model_type, summary) %>%
    unnest
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
# Models ----

#' Run a linear model on terrestrial and epiphytic species
#' separately.
#'
#' @param data Input data. Must include column for "habit"
#' (growth habit, either epiphtyic or terrestrial).
#' @param indep_var Independent variable in data
#' @param resp_var Response variable in data
#'
#' @return Dataframe; model summary
#'
#' @examples
#' # A rather ridiculous example
#' mtcars_hab <- mutate(mtcars, habit = c(
#' rep("terrestrial", 0.5*nrow(mtcars)), 
#' rep("epiphtyic", 0.5*nrow(mtcars)))
#' )
#' run_lm_by_habit(mtcars_hab, "mpg", "wt")
run_lm_by_habit <- function (data, indep_var, resp_var) {
  
  resp_var_enq <- enquo(resp_var)
  indep_var_enq <- enquo(indep_var)
  
  data %>%
    select(y = !!resp_var_enq, x = !!indep_var_enq, habit) %>%
    nest(-habit) %>%
    mutate(
      resp_var = resp_var,
      model = map(data, ~lm(y ~ x, data = .)),
      summ = map(model, tidy)
    ) %>%
    select(-model, -data) %>%
    unnest %>%
    mutate(term = case_when(
      term == "x" ~ indep_var,
      TRUE ~ term
    ))
}

#' Get linear model fits on terrestrial and epiphytic species
#' separately.
#'
#' @param data Input data. Must include column for "habit"
#' (growth habit, either epiphtyic or terrestrial).
#' @param indep_var Independent variable in data
#' @param resp_var Response variable in data
#'
#' @return Dataframe; model summary
#'
#' @examples
#' # A rather ridiculous example
#' mtcars_hab <- mutate(mtcars, habit = c(
#' rep("terrestrial", 0.5*nrow(mtcars)), 
#' rep("epiphtyic", 0.5*nrow(mtcars)))
#' )
#' fit_lm_by_habit(mtcars_hab, "mpg", "wt")
fit_lm_by_habit <- function (data, indep_var, resp_var) {
  
  resp_var_enq <- enquo(resp_var)
  indep_var_enq <- enquo(indep_var)
  
  data %>%
    select(y = !!resp_var_enq, x = !!indep_var_enq, habit) %>%
    nest(-habit) %>%
    mutate(
      resp_var = resp_var,
      indep_var = indep_var,
      model = map(data, ~lm(y ~ x, data = .)),
      fits = map(model, augment)
    ) %>%
    select(-model, -data) %>%
    unnest %>%
    rename(resp_val = y, indep_val = x)
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

# Plotting ----

make_climate_plot <- function (climate_data, site_data, 
                               fits, summaries, 
                               habit_colors, alpha_val = 0.6) {
  
  climate_site_data <- left_join(climate_data, site_data)
  
  mean_temp_fits <-
    fits %>%
    filter(var == "mean_temp") %>%
    rename(mean_temp = value)
  
  # mean temp subplot
  a <- ggplot(climate_site_data, aes(x = el)) +
    geom_line(data = mean_temp_fits, aes(y = .fitted)) +
    geom_point(aes(y = mean_temp, alpha = is_outlier, color = habit)) +
    scale_color_manual(
      values = habit_colors
    ) +
    scale_alpha_manual(
      values = c("yes" = alpha_val, "no" = 1.0)
    ) +
    labs(
      y = "Mean temp. (°C)",
      x = "Elevation (m)"
    ) +
    standard_theme()
  
  sd_temp_fits <-
    fits %>%
    filter(var == "sd_temp") %>%
    rename(sd_temp = value)
  
  # SD temp subplot
  b <- ggplot(climate_site_data, aes(x = el)) +
    geom_line(data = sd_temp_fits, aes(y = .fitted)) +
    geom_point(aes(y = sd_temp, alpha = is_outlier, color = habit)) +
    scale_color_manual(
      values = habit_colors
    ) +
    scale_alpha_manual(
      values = c("yes" = alpha_val, "no" = 1.0)
    ) +
    labs(
      y = "SD Temp. (°C)",
      x = "Elevation (m)"
    ) +
    standard_theme()
  
  min_RH_fits <-
    fits %>%
    filter(var == "min_RH") %>%
    rename(min_RH = value)
  
  # min RH subplot
  c <- ggplot(climate_site_data, aes(x = el, color = habit)) +
    geom_line(data = min_RH_fits, aes(y = .fitted)) +
    geom_point(aes(y = min_RH, alpha = is_outlier)) +
    scale_color_manual(
      values = habit_colors
    ) +
    scale_alpha_manual(
      values = c("yes" = alpha_val, "no" = 1.0)
    ) +
    labs(
      y = "Min. RH (%)",
      x = "Elevation (m)"
    ) +
    standard_theme()
  
  sd_RH_fits <-
    fits %>%
    filter(var == "sd_RH") %>%
    rename(sd_RH = value)
  
  # SD RH subplot
  d <- ggplot(climate_site_data, aes(x = el, color = habit)) +
    geom_line(data = sd_RH_fits, aes(y = .fitted)) +
    geom_point(aes(y = sd_RH, alpha = is_outlier)) +
    scale_color_manual(
      values = habit_colors
    ) +
    scale_alpha_manual(
      values = c("yes" = alpha_val, "no" = 1.0)
    ) +
    labs(
      y = "SD RH (%)",
      x = "Elevation (m)"
    ) +
    standard_theme()
  
  a <- a + blank_x_theme() 
  b <- b + blank_x_theme() 
  
  a + b + c + d + plot_layout(ncol = 2, nrow = 2) &
    theme(legend.position = "none")
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


#' Make scatterplot for two variables by growth
#' habit
#'
#' @param model_fits Results of fitting linear models;
#' output of fit_lm_by_habit
#' @param model_summaries Model summaries; output of
#' run_lm_by_habit
#' @param x_var Independent variable, e.g. "el"
#' @param y_var Dependent variable, e.g. "ntaxa"
#' @param habit_colors Named vector; colors to use
#' to distinguish growth habit
#'
#' @return ggplot object
make_scatterplot_by_habit <- function (model_fits, model_summaries, 
                                       x_var, y_var,
                                       habit_colors) {
  
  model_summaries <-
    model_summaries %>%
    filter(term != "(Intercept)") %>%
    rename(indep_var = term)
  
  plot_data <- 
    left_join(
      model_fits, model_summaries 
    ) %>% filter(resp_var == y_var, indep_var == x_var)
  
  ggplot(plot_data, aes(x = indep_val, color = habit)) +
    geom_line(data = filter(plot_data, p.value < 0.05), aes(y = .fitted)) +
    geom_point(aes(y = resp_val)) +
    labs(
      y = y_var,
      x = x_var
    ) +
    scale_color_manual(
      values = habit_colors
    ) +
    jntools::standard_theme() +
    theme(legend.position = "none")
}
