# Set up captions

# - Figures
fig_nums <- captioner::captioner(prefix = "Fig.")
fig_nums(name = "traits-used-image", caption = "Examples of traits used in this study.")
fig_nums(name = "climate", caption = "Microclimatic variables.")
fig_nums(name = "pca", caption = "Trait PCA.")
fig_nums(name = "tree", caption = "Traits on tree.")
fig_nums(name = "heatmap", caption = "Variable importance scores.")
fig_nums(name = "cwm-div", caption = "Community-weighted means.")
fig_nums(name = "comm-div", caption = "Community diversity metrics.")

# - Tables
table_nums <- captioner::captioner(prefix = "Table")
table_nums(name = "traits-used", caption = "Fern traits used in this study")
table_nums(name = "phylosig-cont", caption = "Phylogenetic signal in continuous traits of ferns from Moorea, French Polynesia")
table_nums(name = "phylosig-binary", caption = "Phylogenetic signal in binary traits of ferns from Moorea, French Polynesia")
table_nums(name = "pics", caption = "Results of phylogenetically independent contrasts analysis of quantitative (sporophyte) traits related to epiphytic growth in ferns from Moorea, French Polynesia")
table_nums(name = "corr-evo", caption = "Pagel's (1994) test of correlated evolution between binary traits in ferns from Moorea, French Polynesia")

# - SI figures
s_fig_nums <- captioner::captioner(prefix = "Fig. S", auto_space = FALSE, suffix = ": ")
s_fig_nums(name = "cwsd", "Comparison of community-weighted standard deviations (CWSD) in traits of epiphytic and terrestrial fern communities on Moorea, French Polynesia")
s_fig_nums(name = "corr-plot", "Correlation plot comparing functional diversity of fern communities on Moorea, French Polynesia calculated using different sets of traits as input")

# - SI tables
s_table_nums <- captioner::captioner(prefix = "Table S", auto_space = FALSE, suffix = ": ")
s_table_nums(name = "species-list", "List of ferns from Moorea, French Polynesia included in this study")
s_table_nums(name = "trait-values", "Trait values of ferns on Moorea, French Polynesia")
s_table_nums(name = "ancova", "Analysis of covariance (ANCOVA) for differences between climatic variables measured with dataloggers placed on the ground (terrestrial growth habit) vs. mounted at 2 m on trees (epiphytic growth habit) on Moorea, French Polynesia, with elevation as a covariate")
s_table_nums(name = "phylosig-binary-strict", "Phylogenetic signal in binary traits of ferns on Moorea, French Polynesia, strict dataset excluding any species whose traits were scored based on taxonomy")
s_table_nums(name = "corr-evo-strict", "Pagel's (1994) test of correlated evolution between binary traits of ferns on Moorea, French Polynesia, strict dataset excluding any species whose traits were scored based on taxonomy")
s_table_nums(name = "model-fit", "Best-fitting models explaining community diversity metrics of ferns on Moorea, French Polynesia by climate and growth habit")

# Make short versions of citation functions
# - Just the number
figure <- pryr::partial(fig_nums, display = "cite")
table <- pryr::partial(table_nums, display = "cite")
s_figure <- pryr::partial(s_fig_nums, display = "cite")
s_table <- pryr::partial(s_table_nums, display = "cite")
# - Just the caption
figure_cap <- function(x) {fig_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
table_cap <- function(x) {table_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
s_figure_cap <- function(x) {s_fig_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
s_table_cap <- function(x) {s_table_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
