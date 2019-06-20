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
  
  get_label <- function(axis_labels, analysis_type_select, axis_select) {
    axis_labels %>% 
      filter(analysis_type == analysis_type_select, axis == axis_select) %>% 
      pull(label)
  }
  
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
