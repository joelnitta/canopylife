plan <- drake_plan(
  
  # Load data ----
  
  # Trait data
  # - Raw specific leaf area (SLA) measurments,
  # includes multiple measures per specimen
  sla_raw = mooreaferns::sla.raw %>% as_tibble,
  # - Other raw trait data, includes one measure per specimen
  morph_raw = mooreaferns::morph.raw %>% as_tibble,
  # - Pre-processed data
  fern_traits = mooreaferns::fern_traits %>% as_tibble,
  
  # Site data
  moorea_sites = mooreaferns::sites %>%
    as_tibble %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),

  # Climate data
  # - Raw data (rel. humidity and temperature every 15 min.)
  moorea_climate_raw = mooreaferns::moorea_climate %>% 
    as_tibble %>% 
    reformat_raw_climate,
  
  # - Process climate data
  daily_mean_climate = get_daily_means(moorea_climate_raw),
  grand_mean_climate = get_grand_means(daily_mean_climate),
  
  # Phylogeny
  phy = mooreaferns::fern_tree,
  
  # Format trait data for SI.
  # Include SD and n for all quantitative (continuous) data
  trait_data_for_si = process_trait_data_for_si(
    sla_raw = sla_raw,
    morph_raw = morph_raw,
    fern_traits = fern_traits
  ),
  
  # Climate analysis ----
  
  # Make full set of climate models including
  # each climate var by elevation, growth habit, and
  # their interaction, then choose the best model for each.
  climate_models = make_climate_models(
    climate_data = grand_mean_climate, 
    site_data = moorea_sites
  ),
  
  # Extract fits of supported models.
  climate_model_fits = extract_model_fits(
    supported_models = climate_models
  ),
  
  # Get summaries of supported models.
  climate_model_summaries = extract_model_summaries(
    supported_models = climate_models
  ),
  
  # Principal components analysis ----
  pca_results = run_trait_PCA(
    traits = fern_traits,
    phy = phy
  ),
  
  # Plots ----
  
  # Set color scheme
  # Epiphytes in green, terrestrial in brown
  habit_colors = brewer.pal(9, "Set1")[c(3,7)] %>%
    set_names(levels(moorea_climate_raw$habit)),
  
  # Make climate plot
  climate_plot = make_climate_plot(
    climate_data = grand_mean_climate,
    site_data = moorea_sites,
    fits = climate_model_fits, 
    summaries = climate_model_summaries, 
    habit_colors = habit_colors)
  
)