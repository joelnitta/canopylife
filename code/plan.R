# Setup workflow plan
plan <- drake_plan(
  
  # Load data ----
  
  # Trait data
  # - Raw specific leaf area (SLA) measurments,
  # includes multiple measures per specimen
  sla_raw = mooreaferns::sla.raw %>% as_tibble,
  
  # - Other raw trait data, includes one measure per specimen
  morph_raw = mooreaferns::morph.raw %>% as_tibble,
  
  # - Pre-processed trait data
  fern_traits = mooreaferns::fern_traits %>% 
    as_tibble %>%
    # make sure binary traits are coded as character
    mutate_at(vars(glands, hairs, gemmae), as.character),
  
  # Site data
  moorea_sites = mooreaferns::sites %>%
    as_tibble %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),
  
  # Climate data
  # - Raw climate data (rel. humidity and temperature every 15 min.)
  moorea_climate_raw = mooreaferns::moorea_climate %>% 
    as_tibble %>% 
    reformat_raw_climate,
  
  # - Process climate data
  daily_mean_climate = get_daily_means(moorea_climate_raw),
  grand_mean_climate = get_grand_means(daily_mean_climate),
  
  # Phylogeny
  phy = mooreaferns::fern_tree,
  
  # Community data for ferns of Moorea only 
  # with rows as species and sites as columns
  comm = mooreaferns::sporocomm %>%
    rownames_to_column("site") %>%
    as_tibble %>%
    filter(site %in% moorea_sites$site) %>%
    gather(species, abundance, -site) %>%
    spread(site, abundance),
  
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
  
  # Phylogenetic signal ----
  
  # - Continuous traits (Blomberg's K and Pagel's lambda)
  phylosig_cont_traits = list(
    as.list(c("stipe", "length", "width", "dissection", "pinna", "sla", "rhizome")),
    traits = list(fern_traits),
    phy = list(phy)
  ) %>%
    pmap_dfr(analyze_cont_phylosig),
  
  # - Binary traits (Fritz and Purvis' D)
  phylosig_binary_traits = analyze_binary_phylosig(
    traits = fern_traits,
    phy = phy
  ),
  
  # Correlated evolution ----
  correlated_evo_test = list(
    as.list(c("morphotype", "glands", "hairs", "gemmae")),
    traits = list(fern_traits),
    phy = list(phy)
  ) %>%
    pmap_dfr(test_corr_evo),
  
  # Phylogenetically independent contrasts ----
  pic_results = run_pic(
    traits = fern_traits,
    phy = phy
  ),
  
  # Phylogenetic community diversity ----
  comm_struc = analyze_comm_struc_by_habit(
    comm = comm, 
    phy = phy, 
    traits = fern_traits
  ),
  
  # Trait community diversity ----
  func_div_epi = analyze_fd_by_habit(
    traits = fern_traits,
    comm = comm,
    habit_type = "epiphytic"
  ),
  
  func_div_ter = analyze_fd_by_habit(
    traits = fern_traits,
    comm = comm,
    habit_type = "terrestrial"
  ),
  
  func_div = bind_rows(func_div_epi, func_div_ter),
  
  # Models ----
  
  # Merge all diversity metrics.
  diversity_data = list(comm_struc, func_div, grand_mean_climate, moorea_sites) %>% 
    reduce(left_join) %>%
    mutate(habit = fct_relevel(habit, c("terrestrial", "epiphytic"))),
  
  # Make list of arguments (dep. var. and indep. var.) for linear models.
  args = cross_df(list(
    resp_vars = c(
      "ntaxa", "mpd.obs.z", "mntd.obs.z", "FDiv", "FEve", "FRic",
      "stipe", "length", "width", "dissection", "pinna", "sla", "rhizome"), 
    indep_vars = c("el", "min_RH")
  )),
  
  # Run linear models on each combination of variables for
  # terrestrial and epiphytic communities separately.
  model_summaries = map2_df(
    .x = args$indep_vars, 
    .y = args$resp_vars, 
    ~ run_lm_by_habit(diversity_data, indep_var = .x, resp_var = .y)
  ),
  
  model_fits = map2_df(
    .x = args$indep_vars, 
    .y = args$resp_vars, 
    ~ fit_lm_by_habit(diversity_data, indep_var = .x, resp_var = .y)
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
    habit_colors = habit_colors),
  
  # Make PCA plot
  pca_plot = make_pca_plot(
    pca_results = pca_results, 
    habit_colors = habit_colors, 
    traits = fern_traits
  )
  
)