This canopylife README.txt file was generated on 2020-01-08 by Joel Nitta

------------------- GENERAL INFORMATION -----------------

Title of Dataset: Life in the canopy: community trait assessments reveal
substantial functional diversity among fern epiphytes

Author Information

	Principal Investigator: Joel H. Nitta Department of Botany, National Museum
	of Natural History, Smithsonian Institution, Washington, D.C., 20013,
	U.S.A. joelnitta@gmail.com

	Associate or Co-investigator: James E. Watkins, Jr. Department of Biology,
	Colgate University, 13 Oak Drive, Hamilton, New York 13346, U.S.A.
	jwatkins@colgate.edu

	Associate or Co-investigator: Charles Davis Department of Organismic and
	Evolutionary Biology and Harvard University Herbaria, Harvard University,
	Cambridge, MA 02138, U.S.A. cdavis@oeb.harvard.edu

Date of data collection: 2012-2015

Geographic location of data collection: Moorea, French Polynesia

Information about funding sources or sponsorship that supported the
collection of the data: Staff at the University of California, Berkeley
Richard B. Gump South Pacific Research Station, Moorea, French Polynesia
provided support with obtaining research and collection permits, laboratory
space, and logistics. Funding provided in part by National Science Foundation
(Doctoral Dissertation Improvement Grant DEB-1311169 to J.H.N. and C.C.D.),
Setup Funds from Harvard University to C.C.D., the American Society of Plant
Taxonomists (Research Grant for Graduate Students to J.H.N.), the Garden Club
of America (Award in Tropical Botany to J.H.N.), the Harvard University
Herbaria (Fernald Fieldwork Fellowship to J.H.N.), the Society of Systematic
Biologists (Graduate Student Research Award to J.H.N.), and the Systematics
Association (Systematics Research Fund to J.H.N.).

--------------------------
SHARING/ACCESS INFORMATION
--------------------------

Licenses/restrictions placed on the data, or limitations of reuse: CC0 1.0
Universal (CC0 1.0)

Recommended citation for the data: Nitta JH, Watkins JEW, Davis CC (2020)
Life in the canopy: community trait assessments reveal substantial functional
diversity among fern epiphytes. Dryad Digital Repository. ADD_DOI_WHEN_AVAILABLE

Citation for and links to publications that cite or use the data: Nitta JH,
Watkins JEW, Davis CC (2020) Life in the canopy: community trait assessments
reveal substantial functional diversity among fern epiphytes. New
Phytologist. ADD_DOI_WHEN_AVAILABLE

Code for analyzing the data is available on github:
https://github.com/joelnitta/canopylife

--------------------
DATA & FILE OVERVIEW
--------------------

File list (filenames, directory structure (for zipped files) and brief
description of all data files):

filmy_SLA.csv: Mass of leaf fragments of (mostly) filmy ferns (family
Hymenophyllaceae) from Moorea, French Polynesia used to calculate specific
leaf area (SLA).

hobo_moorea_aorai_2-24-15.csv: Temperature (celsius) and relative humidity
(%) measured every 15 minutes (Moorea) or 60 minutes (Tahiti) with 33
dataloggers on Moorea and Tahiti, French Polynesia from 2012-07-18 to
2015-02-06.

moorea_climate.csv: Temperature (celsius) and relative humidity (%) measured
every 15 minutes, subsetted to only Moorea and excluding days missing data.

morph_measurements.csv: Measurements of stipe length, frond length, frond
width, rhizome diameter (all in cm), and number of pinna pairs on ferns (one
measurement per individual) from Moorea, French Polynesia.

morph_qual_traits.csv: Observations of qualitative trait states in ferns from
Moorea, French Polynesia.

ppgi_taxonomy.csv: Pteridophyte Phylogeny I working group taxonomic system
for pteridophytes at the genus level and above.

SLA_measurements.csv: Mass of leaf punches of ferns from Moorea, French
Polynesia used to calculate specific leaf area (SLA).

table_1.csv: Summary of traits used in this study.

Additional related data collected that was not included in the current data
package: The following data from Nitta et al. 2017 are already in Dryad
(https://doi.org/10.5061/dryad.df59g) and are used in the analysis:

  - all_plots.csv: Community data matrix including fern sporophytes and
  gametophytes for ferns on Moorea and Tahiti from Nitta et al. 2017
  (data_and_scripts/Comm_Phylo_Analysis/data/all_plots.csv). Rows are sites;
  columns are species. _G after a row name indicates gametophytes; _S
  indicates sporophytes. For sporophytes, values in cells are an abundance
  ranking from 0 to 24 or 25 indicating number of subplots where that species
  was present at a particular site: 0 is completely absent from the site, and
  24 or 25 is maximal abundance (present in all subplots). For gametophytes,
  values in cells are number of individuals observed per species at a
  particular site. Gametophytes identified by BLAST query.

  - treepl_Moorea_Tahiti.tre: Ultramteric, maximum-likelihood tree of ferns on
  Moorea and Tahiti (Newick format) from Nitta et al. 2017
  (data_and_scripts/Comm_Phylo_Analysis/data/treepl_Moorea_Tahiti.tre).

  - species.csv: Taxonomy of all fern and lycophyte species on Moorea and Tahiti
  from Nitta et al. 2017 (data_and_scripts/shared_data/species.csv). "Tahiti
  Only" column (1/0) indicates that a species occurs only on Tahiti and not
  Moorea. "Include" column (1/0) indicates if the species should be included
  in this study or not.

  - sites.csv: Elevation and GPS coordinates (WGS84 datum, decimal format) for each
  site from Nitta et al. 2017 (data_and_scripts/shared_data/sites.csv).
  Columns as follows, site: site name, lat: latitude, long: longitude, el:
  elevation (m).

--------------------------
METHODOLOGICAL INFORMATION
--------------------------

Description of methods used for collection/generation of data: Morphological
traits of ferns in Moorea, French Polynesia were measured from voucher
specimens. Environmental data (relative humidity and temperature) were
measured with data loggers. For additional details, see Nitta, Watkins, and
Davis (2020) New Phyt. ADD_DOI_WHEN_AVAILABLE.

--------------------------
DATA-SPECIFIC INFORMATION
--------------------------

filmy_SLA.csv: Mass of leaf fragments of (mostly) filmy ferns (family
Hymenophyllaceae) from Moorea, French Polynesia used to calculate specific
leaf area (SLA). These specimens had leaves with veins that were too close to
avoid when preparing leaf punches (see SLA_measurements.csv), so irregularly
shaped leaf fragments avoiding veins were cut out from dried leaf laminae and
used instead. Area measured with ImageJ (Abramoff 2014) using digital images
of leaf fragments.

Number of variables: 9

Number of cases/rows: 155

Variable list, defining any abbreviations, units of measure, codes or symbols
used:
	- species: Species name.
	- specimen: Voucher number (J.H. Nitta collection number).
	- individual: Arbitrarily assigned number for each individual per specimen.
	- mass_g: Mass of the leaf fragment in (g).
	- area_cm2: Area of the leaf fragment (square cm).
	- sla_m2_per_kg: Specific leaf area (square m per kg).
	- notes: Miscellaneous notes.
	- comments: Additional miscellaneous notes.
	- original_weight: Original weight if leaf fragment was re-weighed.

Missing data codes: Missing data have no values (nothing entered between
commas in the CSV file).

Specialized formats or other abbreviations used: None.

--------------------------

hobo_moorea_aorai_2-24-15.csv: Temperature (celsius) and relative humidity
(%) measured every 15 minutes (Moorea) or 60 minutes (Tahiti) with 33
dataloggers on Moorea and Tahiti, French Polynesia from 2012-07-18 to
2015-02-06. On Moorea, Hobo Pro v2 data loggers with the RS3 Solar Radiation
Shield (Onset Corporation, Bourne, Massachusetts, USA) were used. On Tahiti,
RHTemp 1000 data loggers (MadgeTech, Warner, New Hampshire, USA) protected
with custom radiation shields made from plastic circuit boxes were used.

Number of variables: 5

Number of cases/rows: 1464604

Variable list:
	- date: Date (YYYY-MM-DD).
	- time: Time (hh::mm::ss, 24 hr format).
	- temp: Temperature (celsius).
	- RH: Relative humidity (percent).
	- site: Site name and placement of datalogger. Site name formatted as name
	  of mountain (e.g., "Rotui") then approximate elevation (e.g., "400m")
	  separated by underscore. Placement of datalogger given by "ter" (at
	  ground level) or "epi" (mounted ca. 2m height on a tree) after the site
	  name separated by an underscore. GPS positions of sites can be found in
	  Nitta et al. 2017.

Missing data codes: No missing data.

Specialized formats or other abbreviations used: None.

--------------------------

moorea_climate.csv: Microclimate data (same data as
hobo_moorea_aorai_2-24-15.csv), subsetted to only Moorea and excluding days
missing data. Dates range from 2013-07-07 to 2014-07-05. The terrestrial
datalogger at the "tohiea_400m" site malfunctioned and was missing data for a
substantial part of the final survey period (2014-03-12 to 2014-07-05), so
data for the same dates from 2013 were recoded as 2014 and used instead for
this datalogger only.

Number of variables: 7

Number of cases/rows: 587,912

Variable list:
  - date: Date (YYYY-MM-DD).
  - time: Time (hh::mm::ss, 24 hr format).
  - site: Site name, formatted as name of mountain (e.g., "Rotui") then
    approximate elevation (e.g., "400m") separated by underscore. GPS
    positions of sites can be found in Nitta et al. 2017.
  - habit: Placement of datalogger given by "ter" (at ground level) or "epi"
    (mounted ca. 2m height on a tree)
	- temp: Temperature (celsius).
	- rh: Relative humidity (percent).
  - vpd: Vapor-pressure deficit (kPa). Calculated from relative humidity and
    temperature using plantecophys::RHtoVPD().

Missing data codes: No missing data.

Specialized formats or other abbreviations used: None.

--------------------------

morph_measurements.csv: Measurements of stipe length, frond length, frond
width, rhizome diameter (all in cm), and number of pinna pairs on ferns (one
measurement per individual) from Moorea, French Polynesia.

Number of variables: 11

Number of cases/rows: 261

Variable list:
  - specimen: Voucher number (J.H. Nitta collection number). If specimen is
	"1", that means no voucher specimen was prepared when the measurement was
	made.
  - species: Species name.
  - stipe_length: Length of stipe (petiole) (cm).
  - frond_length: Length of frond (leaf) (cm).
  - lamina_width: Width of lamina at widest point (cm).
  - rhizome_dia: Diameter of rhizome at widest point (cm).
  - pinna_pairs: Number of pinnae pairs.
  - source: Source of data. "measurement" indicates data measured from
	specimen; other values indicate name of source of data. For full references
	see, Supporting Information in Nitta et al. (2020).
  - comments: Miscellaneous notes.
  - notes: Additional miscellaneous notes.
  - exclude: Should the measurement be excluded from analysis? "yes" if true,
	otherwise missing.

Missing data codes: Missing data have no values (nothing entered between
commas in the CSV file).

Specialized formats or other abbreviations used: None.

--------------------------

morph_qual_traits.csv: Observations of qualitative trait states in ferns from
Moorea, French Polynesia.

Number of variables: 9

Number of cases/rows: 143

Variable list:
  - species: Species name.
  - dissection: Integer from one to 10 describing the degree of lamina
    dissection as follows: 1 = simple, 2 = pinnatifid or pinnatisect, 3 =
    1-pinnate, 4 = 1-pinnate-pinnatifid, 5 = 2-pinnate, 6 =
    2-pinnate-pinnatifid, 7 = 3-pinnate, 8 = 3-pinnate-pinnatifid, 9 = more
    than 3-pinnate-pinnatifid, 10 = binpinnatifid or tripinnatifid.
  - habit_binary: Growth habit coded as terrestrial (0) or epiphytic (1);
	  epiphytes defined by lacking a connection to the soil.
  - glands: Presence (1) or absence (0) of glands on gametophytes.
  - hairs: Presence (1) or absence (0) of hairs on gametophytes.
  - morphotype: Gametophyte morphotype sensu Farrar, et al. 2008
  - gemmae: Presence (1) or absence (0) of gemmae on gametophytes.
  - source: Data source for gametophyte traits, coded as follows: 1 = lab
	  observation, 2 = Nayar and Kauar 1971, 3 = Lloyd 1980, 4 = field
	  observation, 5 = genus/family level character, 6 = Zhang et al 2008, 7 =
	  Bierhorst 1967, 8 = Martin et al 2006, 9 = Tigerschiöld 1989, 10 =
	  Tigerschiöld 1990, 11 = Atkinson 1975, 12 = Chen 2014. For full
		references, see Supporting Information in Nitta et al. (2020).
  - notes: Miscellaneous notes.

Missing data codes: Missing data have no values (nothing entered between
commas in the CSV file).

Specialized formats or other abbreviations used: None.

--------------------------

ppgi_taxonomy.csv: Pteridophyte Phylogeny I working group taxonomic system
for pteridophytes at the genus level and above (Pteridophyte Phylogeny Group
I, 2016).

Number of variables: 8

Number of cases/rows: 337

Variable list:
	- class: Name of class.
	- order_id: Unique ID for order.
	- order: Name of order.
	- suborder: Name of suborder.
	- fam_id: Unique ID for family.
	- family: Name of family.
	- subfamily: Name of subfamily.
	- genus: Name of genus.

Missing data codes: Non-applicable values coded as "NA".

Specialized formats or other abbreviations used: none.

--------------------------

SLA_measurements.csv: Mass of leaf punches of ferns from Moorea, French
Polynesia used to calculate specific leaf area (SLA). Punches made with
biopsy punches on fresh leaves then dried overnight in a lab heat oven. Leaf
laminae excluding veins was used as much as possible; 2 mm punches were used
instead of 4 mm when needed to avoid veins.

Number of variables: 7

Number of cases/rows: 1914

Variable list:
  - order: Unique row number (1-1914).
  - specimen: Specimen ID number. Four-digit numbers are J.H. Nitta
    collection numbers (voucher specimen numbers). Values formatted as "SLA
    (year) sample (number)" indicate leaf samples without any corresponding
    voucher specimen. Values with "SKIP" are row placeholders that don't
    contain any data.
  - species: Species name.
  - mass: Mass of leaf punch (g). If >1 punch was taken (as indicated in
    number_of_punches), this is the total mass of all punches for that
    sample.
  - notes: Miscellaneous notes.
  - number_of_punches: Number of punches. If empty, this indicates 1 punch.
  - size_of_punches: Diameter of punches (mm).

Missing data codes: Missing data have no values (nothing entered between
commas in the CSV file). Note that for number_of_punches, an empty value
indicates 1 punch.

Specialized formats or other abbreviations used: None.

--------------------------

table_1.csv: Summary of traits used in this study. Used to produce Table 1 in
Nitta et al. 2020.

Number of variables: 4

Number of cases/rows: 11

Variable list:
  - trait: Name of trait.
  - data_type: Type of data and unit, if applicable.
  - ecological_significance: Ecological significance.
  - reference: Reference for ecological significance.

Missing data codes: No missing data.

Specialized formats or other abbreviations used: None.

--------------------------
REFERENCES
--------------------------

Abràmoff, M. D., P. J. Magalhães, and S. J. Ram. 2004. Image processing with
imageJ. Biophotonics International 11:36–41.

Nitta, J. H., J.-Y. Meyer, R. Taputuarai, and C. C. Davis. 2017. Life cycle
matters: DNA barcoding reveals contrasting community structure between fern
sporophytes and gametophytes. Ecological Monographs 87:278–296.

Nitta, J. H., J. E. Watkins, and C. C. Davis. 2020. Life in the canopy:
community trait assessments reveal substantial functional diversity among
fern epiphytes. New Phytologist. ADD_DOI_WHEN_AVAILABLE

Pteridophyte Phylogeny Group I. 2016. A community-derived classification for
extant lycophytes and ferns. Journal of Systematics and Evolution 54:563–603.
