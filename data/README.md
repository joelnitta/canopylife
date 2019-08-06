# Data

Pre-processed data files.

table_1.csv: Table of traits used in this study.

ppgi_taxonomy.csv: Taxonomy of ferns and lycophytes following Pteridophyte Phylogeny Group I. 2016. A community-derived classification for extant lycophytes and ferns. Journal of Systematics and Evolution 54:563–603.

nitta_2017 contains the following data files from Nitta, et al. 2017. Life cycle matters: DNA barcoding reveals contrasting community structure between fern sporophytes and gametophytes. Ecological Monographs 87:278–296. This folder and its contents will be downloaded and unzipped by running make.R.

- all_plots.csv: Community data matrix including fern sporophytes and gametophytes for ferns on Moorea and Tahiti from Nitta et al. 2017 (data_and_scripts/Comm_Phylo_Analysis/data/all_plots.csv). Rows are sites; columns are species. _G after a row name indicates gametophytes; _S indicates sporophytes. For sporophytes, values in cells are an abundance ranking from 0 to 24 or 25 indicating number of subplots where that species was present at a particular site: 0 is completely absent from the site, and 24 or 25 is maximal abundance (present in all subplots). For gametophytes, values in cells are number of individuals observed per species at a particular site. Gametophytes identified by BLAST query.

- treepl_Moorea_Tahiti.tre: Ultramteric, maximum-likelihood tree of ferns on Moorea and Tahiti (Newick format) from Nitta et al. 2017 (data_and_scripts/Comm_Phylo_Analysis/data/treepl_Moorea_Tahiti.tre).

- species: Taxonomy of all fern and lycophyte species on Moorea and Tahiti from Nitta et al. 2017 (data_and_scripts/shared_data/species.csv). "Tahiti Only" column (1/0) indicates that a species occurs only on Tahiti and not Moorea. "Include" column (1/0) indicates if the species should be included in this study or not.

- sites: Elevation and GPS coordinates (WGS84 datum, decimal format) for each site from Nitta et al. 2017 (data_and_scripts/shared_data/sites.csv). Columns as follows, site: site name, lat: latitude, long: longitude, el: elevation (m). 