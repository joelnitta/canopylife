# Setup workflow plan
plan <- drake_plan(
  
  # Load data ----
  
  ### Trait data ###
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
  
  # - Also make "strict" trait dataset by excluding species with 
  # gameto traits only known from taxonomy
  fern_traits_strict = filter(fern_traits, gameto_source != "T"),
  
  ### Site data with elevation ###
  moorea_sites = mooreaferns::sites %>%
    as_tibble %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),
  
  ### Climate data ###
  # - Raw climate data (rel. humidity and temperature every 15 min.)
  # - Also calculate vapor pressure deficit (VPD) from temp and RH
  moorea_climate_raw = mooreaferns::moorea_climate %>% 
    as_tibble %>% 
    reformat_raw_climate,
  
  # - Process climate data
  daily_mean_climate = get_daily_means(moorea_climate_raw),
  grand_mean_climate = get_grand_means(daily_mean_climate),
  
  ### Phylogeny ###
  phy = mooreaferns::fern_tree,
  
  ### Community data for ferns of Moorea ###
  # (with rows as species and sites as columns)
  comm = mooreaferns::sporocomm %>%
    rownames_to_column("site") %>%
    as_tibble %>%
    filter(site %in% moorea_sites$site) %>%
    gather(species, abundance, -site) %>%
    spread(site, abundance),
  
  ### PPGI taxonomy ###
  ppgi = read_csv(file_in("data/ppgi_taxonomy.csv")),
  
  ### Format trait data for SI ###
  # Include SD and n for all quantitative (continuous) data
  trait_data_for_si = process_trait_data_for_si(
    sla_raw = sla_raw,
    morph_raw = morph_raw,
    fern_traits = fern_traits
  ),
  
  # Climate analysis ----
  
  ### Prepare data ###
  # Subset climate variables to only those with correlation
  # coefficients less than 0.9,
  # and add elevation.
  climate_select = select_climate_vars(grand_mean_climate) %>%
    left_join(moorea_sites),
  
  # Make three versions of climate data for further analysis:
  # - Square-root transform selected variables
  # so they have a more normal distribution.
  climate_trans = climate_select %>%
    mutate_at(vars(temp_sd, vpd_mean, vpd_min, vpd_sd), sqrt),
  
  # - Remove the outliers, but don't transform.
  climate_no_outlier = climate_select %>%
    filter(is_outlier == "no"),
  
  # - Remove outliers and transform.
  climate_trans_no_outlier = climate_trans %>%
    filter(is_outlier == "no"),
  
  # - Make vector of climate variables for modeling
  climate_vars = climate_select %>%
    select(-site, -habit, -el, -is_outlier, -lat, -long) %>%
    colnames() %>%
    rlang::set_names(),
  
  ### Run ANCOVA ###
  climate_ancova_results = run_ancova(
    data = climate_trans_no_outlier, 
    resp_vars = climate_vars),
  
  ### Make linear models ###
  
  # - Make full set of climate models including
  # each climate var by elevation, growth habit, and
  # their interaction, then choose the best model for each.
  # Exclude single outlier site, but use untransformed
  # data so the fits can be plotted.
  climate_models = choose_habit_elevation_models(
    data = climate_trans_no_outlier,
    resp_vars = climate_vars),
  
  # - Extract fits of best models.
  # Transform fit of square-rooted variables back for plotting.
  climate_model_fits = extract_model_fits(climate_models) %>%
    mutate(.fitted = case_when(
      var %in% c("temp_sd", "vpd_mean", "vpd_min", "vpd_sd") ~ .fitted^2,
      TRUE ~ .fitted)
    ),
  
  # - Get summaries of best model parameters.
  climate_model_parameters = extract_model_parameters(climate_models),
  
  # - Get summaries of best models.
  climate_model_summaries = extract_model_summaries(climate_models),
  
  # - Check for spatial autocorrelation in residuals.
  climate_model_moran = climate_models %>%
    select(-model) %>%
    unnest %>%
    select(var, model_type, AICc, morans_I = statistic, I_pval = p.value),
  
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
  
  # - also run for "strict" trait dataset
  phylosig_binary_traits_strict = analyze_binary_phylosig(
    traits = fern_traits_strict,
    phy = phy
  ),
  
  # Correlated evolution ----
  correlated_evo_test = list(
    as.list(c("morphotype", "glands", "hairs", "gemmae")),
    traits = list(fern_traits),
    phy = list(phy)
  ) %>%
    pmap_dfr(test_corr_evo),
  
  # - also run for "strict" trait dataset
  correlated_evo_test_strict = list(
    as.list(c("morphotype", "glands", "hairs", "gemmae")),
    traits = list(fern_traits_strict),
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
  # (includes quantitative, i.e., sporophyte, traits only).
  # Run separately for epiphytes and terrestrial then combine.
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
  
  # Modeling community diversity ----
  
  ### Prepare diversity data (response variables) ###
  
  # - Combine diversity metrics into single dataframe
  div_metrics_all = left_join(comm_struc, func_div) %>%
    select(site, habit,
           ntaxa, mpd.obs.z, mntd.obs.z, 
           FRic, FEve, FDiv,
           dissection, sla, stipe, length, width, rhizome, pinna),
  
  # - Subset diversity metrics to only those with correlation
  # coefficients less than 0.9, and add elevation and outlier status.
  div_metrics_select = select_div_metrics(div_metrics_all) %>%
    left_join(moorea_sites) %>%
    mutate(
      is_outlier = case_when(
        site == "Rotui_800m_slope" ~ "yes",
        TRUE ~ "no"
      )
    ),
  
  ### Run full-subset analysis using GAMs ###
  
  # - Merge all diversity metrics and climate data for GAMs
  # (transformed version of climate, but including outlier)
  div_climate_trans = inner_join(div_metrics_select, climate_trans),
  
  # - Make named vector of response variables to test
  resp_vars = div_metrics_select %>%
    select(-site, -habit, -el, -is_outlier, -long, -lat) %>% 
    colnames %>% rlang::set_names(),
  
  # - Run full-subsets model analysis.
  fss_div_results = purrr::map(
    resp_vars, 
    ~ run_full_subset_canopy_mods(
      data = div_climate_trans,
      indep_vars = climate_vars,
      resp_var = .)
    ),
  
  # - Extract table of importance of each environmental variable.
  important_div_vars = get_important_vars(fss_div_results),
  
  # - Extract table of best-fit models.
  best_fit_div_models = get_best_fss_mods(fss_div_results),
  
  # - Test for spatial autocorrelation in residuals with Moran's I.
  best_fit_div_models_moran = mutate(
      best_fit_div_models,
      moran_test = map2(
        .x = resp_var, .y = modname,
        ~ check_moran_fss(
          model_set = fss_div_results,
          resp_var = .x,
          best_mod = .y,
          site_data = moorea_sites
        )
      )
    ) %>%
    unnest() %>%
    select(resp_var, modname, AICc, r2, delta_AICc, 
           morans_I = statistic, morans_I_pval = p.value),
  
  ### Run linear models of diversity metrics by elevation ###
  
  # - Make full set of models including each diversity metric 
  # by elevation, growth habit, and
  # their interaction, then choose the best model for each.
  div_el_models = choose_habit_elevation_models(div_metrics_select, resp_vars),
  
  # - Extract fits of best models.
  div_el_model_fits = extract_model_fits(div_el_models),
  
  # - Get summaries of best model parameters.
  div_el_model_parameter_summaries = extract_model_parameters(div_el_models),
  
  # - Get summaries of best models.
  div_el_model_summaries = extract_model_summaries(div_el_models),
  
  # - Check for spatial autocorrelation in residuals.
  div_el_model_moran = div_el_models %>%
    select(-model) %>%
    unnest %>%
    select(var, model_type, AICc, morans_I = statistic, I_pval = p.value),
  
  ### Run t-tests by growth habit ###
  div_t_test_results = map_df(
    resp_vars,
    ~ run_t_test_by_habit(
      resp_var = .,
      data = div_metrics_select
    )
  ),
  
  # Plots ----
  
  # Set color scheme: epiphytes in green, terrestrial in brown.
  habit_colors = brewer.pal(9, "Set1")[c(3,7)] %>%
    set_names(levels(moorea_climate_raw$habit)),
  
  # Make climate plot.
  climate_plot = make_climate_plot(
    data = climate_select,
    resp_vars = climate_vars,
    fits = climate_model_fits, 
    summaries = climate_model_summaries, 
    habit_colors = habit_colors),
  
  # Make PCA plot.
  pca_plot = make_pca_plot(
    pca_results = pca_results, 
    habit_colors = habit_colors, 
    traits = fern_traits
  ),
  
  # Make plot of traits with tree.
  traits_with_tree = plot_traits_on_tree(
    traits = fern_traits,
    phy = phy,
    ppgi = ppgi
  ),
  
  # Make community diversity scatterplots.
  div_scatterplots = map(
    resp_vars, 
    ~ make_scatterplot(
      data = div_metrics_select,
      fits = div_el_model_fits, 
      summaries = div_el_model_summaries, 
      yval = .,
      habit_colors = habit_colors)
  ),
  
  # Make community diversity boxplots.
  div_boxplots = map(
    resp_vars, 
    ~ make_boxplot(
      data = div_metrics_select,
      summaries = div_t_test_results, 
      yval = .,
      habit_colors = habit_colors)
  ),
  
  # Combine community diversity scatterplots
  # and boxplots into final figure.
  combined_comm_div_plots = combine_comm_div_plots(
    scatterplots = div_scatterplots, 
    boxplots = div_boxplots),
  
  # Combine community-weighted means scatterplots
  # and boxplots into final figure.
  combined_cwm_plots = combine_cwm_plots(
    scatterplots = div_scatterplots, 
    boxplots = div_boxplots),
  
  # Make heatmap of importance scores.
  importance_heatmap = make_heatmap(important_div_vars),
  
  # Write out MS ----
  ms = rmarkdown::render(
    knitr_in("ms/manuscript.Rmd"),
    quiet = TRUE),
  
  # Write out supplemental information
  si = rmarkdown::render(
    knitr_in("si/SI.Rmd"),
    quiet = TRUE)
  
)
