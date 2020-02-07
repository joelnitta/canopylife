# Setup workflow plan
plan <- drake_plan(
  
  # Data processing ----
  
  ### Pre-processed data ###
  # Download and unzip data from Nitta et al. 2017 Ecol. Monographs from Dryad
  # The dataset must be downloaded first by going to
  # https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g, clicking on
  # "Download dataset", and saving to the "data" folder in this project.
  nitta_2017_data = unzip_nitta_2017(
    zipped_path = file_in("data/doi_10.5061_dryad.df59g__v1.zip"),
    unzip_path = "data/nitta_2017",
    # Track data files used as input in analyses
    out1 = file_out("data/nitta_2017/sites.csv"),
    out2 = file_out("data/nitta_2017/treepl_Moorea_Tahiti.tre"),
    out3 = file_out("data/nitta_2017/all_plots.csv"),
    out4 = file_out("data/nitta_2017/species.csv")),
  
  ### Trait data ###
  # - Read in list of accepted species to use for this study
  # i.e., ferns of Moorea (129 spp)
  species_list = read_csv(file_in("data/nitta_2017/species.csv")) %>%
    clean_names() %>%
    filter(include == 1, tahiti_only == 0) %>%
    # Exclude Microsorum_xmaximum (hybrid between M. grossum and M. commutatum)
    filter(genus_sp != "Microsorum_xmaximum") %>%
    pull(genus_sp),
  
  # - Process raw specific leaf area (SLA) measurments
  # (includes multiple measures per specimen).
  # SLA is in sq m per kg
  sla_raw = combine_raw_sla(
    file_in("data/SLA_measurements.csv"),
    file_in("data/filmy_SLA.csv"),
    species_list
  ),
  
  # - Process raw continuous trait data
  # (includes one measure per specimen,
  #  multiple specimens per species).
  # All measurements in cm, except number of pinna pairs.
  morph_cont_raw = process_raw_cont_morph(
    file_in("data/morph_measurements.csv"),
    species_list
  ),
  
  # - Process raw qualitative trait data
  # (includes one obsevation per specimen).
  morph_qual_raw = process_raw_qual_morph(
    file_in("data/morph_qual_traits.csv"),
    species_list
  ),
  
  # - Combine raw trait data into final trait matrix
  fern_traits = make_trait_matrix(sla_raw, morph_cont_raw, morph_qual_raw) %>%
    select(-source_1, -source_2),
  
  # - Also make "strict" trait dataset by excluding species with 
  # gameto traits only known from taxonomy
  fern_traits_strict = make_trait_matrix(sla_raw, morph_cont_raw, morph_qual_raw) %>%
    mutate_at(vars(source_1, source_2), ~replace_na(., "none")) %>%
    filter(source_1 != "T", source_2 != "T") %>%
    select(-source_1, -source_2),
  
  # - Also make trait dataset with continuous traits scaled, but not log-transformed.
  fern_traits_scaled = transform_traits(
    fern_traits, 
    log_trans = FALSE, 
    scale_traits = TRUE,
    scale_select = c("sla", "stipe", "length", 
                     "width", "rhizome", "pinna",
                     "dissection")
  ),
  
  # - Also make trait dataset with continuous traits scaled and 
  # size-related traits log-transformed.
  fern_traits_log_scaled = transform_traits(
    fern_traits, 
    log_trans = TRUE, 
    scale_traits = TRUE,
    trans_select = c("length", "rhizome", "stipe", "width"),
    scale_select = c("sla", "stipe", "length", 
                     "width", "rhizome", "pinna",
                     "dissection")
  ),
  
  ### Site data with elevation ###
  moorea_sites = read_csv(
    file_in("data/nitta_2017/sites.csv")) %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),
  
  ### Climate data ###
  
  # - Read in raw climate data (rel. humidity and temperature) from all dataloggers 
  # on Moorea and Tahiti from 2012-07-18 to 2015-02-06.
  climate_raw = read_csv(file_in("data/hobo_moorea_aorai_2-24-15.csv")),
  
  # - Process raw climate data: 
  # subset to only Moorea and days without missing data,
  # replace one chunk of missing data with same dates from previous year,
  # calculate vapor pressure deficit (VPD) from temp and RH.
  moorea_climate = process_raw_climate(climate_raw),
  
  # - Write out moorea_climate as CSV file to include with Dryad data.
  moorea_climate_out = moorea_climate %>%
    select(date, time, site, habit, temp, rh, vpd) %>%
    write_csv(file_out("results/moorea_climate.csv")),
  
  # - Calculate mean climate variables.
  daily_mean_climate = get_daily_means(moorea_climate),
  grand_mean_climate = get_grand_means(daily_mean_climate),
  
  ### Phylogeny ###
  # Ferns on Moorea and Tahiti (146 species)
  phy = ape::read.tree(
    file_in("data/nitta_2017/treepl_Moorea_Tahiti.tre"
    )) %>%
    # Exclude Microsorum_xmaximum (hybrid between M. grossum and M. commutatum)
    drop.tip("Microsorum_xmaximum"),
  
  ### Community data for ferns of Moorea ###
  # Including sporophyte only,
  # with rows as species and sites as columns (130 x 18).
  # (17 sites with 1 column for species name)
  comm = process_community_matrix(
    file_in("data/nitta_2017/all_plots.csv"),
    species_list,
    moorea_sites),
  
  ### PPGI taxonomy ###
  ppgi = read_csv(file_in("data/ppgi_taxonomy.csv")),
  
  # Format data for SI ----
  # - Traits: include SD and n for all quantitative (continuous) data
  trait_data_for_si = process_trait_data_for_si(
    sla_raw, morph_cont_raw, morph_qual_raw
  ),
  
  # - Species names
  # Look up scientific names (specis names with taxonomic author) using
  # a live search of taxonomic databases via the Global Name Resolver 
  # (https://resolver.globalnames.org/). Since results have the possibility
  # of changing, the results are checked-in to the repo at the time of manuscript 
  # submission.
  taxonomic_data_for_si = reformat_species_names(species_list) %>%
    lookup_taxonomy,
  
  # Climate analysis ----
  
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
  # Exclude single outlier site.
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
    traits = fern_traits_log_scaled,
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
  
  # Run Pagel's test of corr. evol. on binary (gametophyte) traits.
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
  
  # For SI, run test on widespread gametophytes vs. pres/abs of gemmae.
  correlated_evo_test_gemmae = test_corr_evo_range_gemmae(
    file_in("data/nitta_2017/all_plots.csv"), 
    phy, 
    fern_traits, 
    moorea_sites
  ),
  
  # Phylogenetically independent contrasts ----
  pic_results = run_pic(
    traits = fern_traits,
    phy = phy
  ),
  
  # Phylogenetic community diversity ----
  phy_comm_struc = analyze_phy_struc_by_habit(
    comm = comm, 
    phy = phy, 
    traits = fern_traits,
    null_model = "independentswap",
    iterations = 1000
  ),
  
  # Functional community diversity ----
  # Define sets of traits for analysis:
  # - all traits
  all = c("stipe", "length", "width", "rhizome",
          "sla", "pinna", "dissection",
          "morphotype", "glands", "hairs", "gemmae"),
  
  # - sporophyte traits related to size
  size = c("stipe", "length", "width", "rhizome"),
  
  # - sporophyte traits other than size
  other = c("sla", "pinna", "dissection"),
  
  # - gametophyte traits
  gameto = c("morphotype", "glands", "hairs", "gemmae"),
  
  # - reduced set of traits excluding frond length and width,
  # which are highly correlated with stipe length.
  reduced = c("stipe", "rhizome",
              "sla", "pinna", "dissection",
              "morphotype", "glands", "hairs", "gemmae"),
  
  # For SI, test out sets of different traits and transformations
  func_comm_struc_test = target(
    analyze_func_struc_by_habit(
      comm = comm,
      traits = trait_data,
      traits_select = traits_select_list,
      null_model = "independentswap",
      abundance_weighted = TRUE,
      iterations = 10000
    ),
    transform = cross(
      trait_data = c(fern_traits_scaled, fern_traits_log_scaled),
      traits_select_list = c(all, size, other, gameto, reduced),
    )
  ),
  
  func_comm_struc_combined = target(
    bind_data(func_comm_struc_test),
    transform = combine(func_comm_struc_test)
  ),
  
  # For main results, analyze traits in the "reduced" set
  # that have been log-transformed and scaled.
  func_comm_struc = analyze_func_struc_by_habit(
    comm = comm,
    traits = fern_traits_log_scaled,
    traits_select = reduced,
    null_model = "independentswap",
    # Weigh so that sporophyte and gametophyte traits
    # are split evenly
    weights = c(rep(.5 / 5, 5), rep(.5 / 4, 4)),
    abundance_weighted = TRUE,
    iterations = 10000
  ),
  
  # Calculate community-weighted means with standard deviation
  # - in long format
  cwm_long = calculate_cwm(fern_traits, comm, moorea_sites),
  # - and by site
  cwm_by_site = cwm_long %>% select(-sd, -el) %>% spread(trait, cwm),
  
  # Modeling ----
  
  ### Prepare diversity data (response variables) ###
  
  # - Combine diversity metrics into single dataframe
  div_metrics_all = left_join(
    select(phy_comm_struc, 
           site, habit, ntaxa, 
           mpd.obs.z.phy = mpd.obs.z, 
           mntd.obs.z.phy = mntd.obs.z),
    select(func_comm_struc,
           site, habit, 
           mpd.obs.z.func = mpd.obs.z, 
           mntd.obs.z.func = mntd.obs.z)) %>%
    left_join(cwm_by_site) %>%
    left_join(moorea_sites),
  
  ### Run full-subset analysis using GAMs ###
  
  # - Merge all diversity metrics and climate data for GAMs
  # (transformed version of climate, but including outlier)
  div_climate_trans = inner_join(div_metrics_all, climate_trans),
  
  # - Make named vector of response variables to test
  resp_vars = div_metrics_all %>%
    select(-site, -habit, -el, -long, -lat) %>%
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
  div_el_models = choose_habit_elevation_models(div_metrics_all, resp_vars),
  
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
      data = div_metrics_all
    )
  ),
  
  # Plots ----
  
  # Set color scheme: epiphytes in green, terrestrial in brown.
  habit_colors = brewer.pal(9, "Set1")[c(3,7)] %>%
    set_names(levels(moorea_climate$habit)),
  
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
  
  # Make CWM scatterplots.
  cwm_scatterplots = map2(
    .x = c("stipe", "length", "rhizome",
           "sla", "pinna", "dissection") %>% set_names(.),
    .y = c("Stipe length (cm)", "Frond length (cm)", "Rhizome dia. (cm)",
           "SLA", "Pinna no.", "Frond dissection"),
    ~ make_cwm_scatterplot(
      data = cwm_long,
      fits = div_el_model_fits,
      summaries = div_el_model_summaries,
      t_test_results = div_t_test_results,
      yval = .x,
      ylab = .y,
      habit_colors = habit_colors)
  ),
  
  # Combine community-weighted means scatterplots into final figure.
  combined_cwm_plots = combine_cwm_plots(
    scatterplots = cwm_scatterplots),
  
  # Make community diversity scatterplots
  div_scatterplots = map(
    c("ntaxa", "mpd.obs.z.phy",  "mntd.obs.z.phy", 
      "mpd.obs.z.func", "mntd.obs.z.func") %>% set_names(.),
    ~ make_div_scatterplot(
      data = div_metrics_all,
      fits = div_el_model_fits,
      summaries = div_el_model_summaries,
      yval = .,
      habit_colors = habit_colors)
  ),
  
  # Make community diversity boxplots.
  div_boxplots = map(
    c("ntaxa", "mpd.obs.z.phy",  "mntd.obs.z.phy",  
      "mpd.obs.z.func", "mntd.obs.z.func") %>% set_names(.),
    ~ make_boxplot(
      data = div_metrics_all,
      summaries = div_t_test_results,
      yval = .,
      habit_colors = habit_colors)
  ),
  
  # Combine community diversity scatterplots
  # and boxplots into final figure.
  combined_comm_div_plots = combine_comm_div_plots(
    scatterplots = div_scatterplots,
    boxplots = div_boxplots),
  
  # Make heatmap of importance scores.
  importance_heatmap = make_heatmap(important_div_vars),
  
  # Make jitter plus box plot of community-wide SD values
  # for SI
  cwsd_scatterplot = make_cwsd_scatterplot(cwm_long, habit_colors),
  
  # Manuscript ----
  
  # First render to PDF, keeping the latex
  ms_pdf = render_tracked(
    knitr_in("ms/manuscript.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results"),
    tracked_output = file_out(here::here("results/manuscript.tex"))
  ),
  
  # Next use the latex to convert to docx with pandoc
  ms_docx = latex2docx(
    latex = file_in(here::here("results/manuscript.tex")),
    docx = file_out(here::here("results/manuscript.docx")),
    template = file_in(here::here("ms/new-phytologist.docx")),
    wd = here::here("results")
  ),
  
  # SI ----
  si_pdf = rmarkdown::render(
    knitr_in("ms/SI.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results")
  )
  
)
