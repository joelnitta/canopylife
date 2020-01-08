# README

Old versions of raw data not used in current project.

## Files

morph_qual_traits_2019-07-29.csv: Most recent qualitative trait data, and the one used for final analysis. Modified from morph_qual_traits_2016-09-03_mod.csv. The exact same as data/morph_qual_traits.

morph_qual_traits_2016-09-03_mod.csv: Older version of qualitative trait modified from 
Moorea_ferns_morph_2016-09-03.xls "morphology" sheet. Includes hard-coded mean 
quantitative values from earlier data (and various other qual traits).

Moorea_ferns_morph_2016-09-03.xls: Final morphology data copied from Moorea CommPhy project.

## Change log

2019-07-29: Add morph_qual_traits_2019-07-29 by copying from
morph_qual_traits_2016-09-03_mod.csv then modifying as follows: 

Delete unused columns, change 'OTU' to 'species'.

Add the following sources:
- "Acrophorus_raiateensis" <- "4"
- "Antrophyum_sp1" <- "4"
- "Ophioglossum_pendulum" <- "2"
- "Oreogrammitis_subspathulata" <- "4"
- "Schizaea_fistulosa" <- "2"
- "Scleroglossum_sulcatum" <- "5"
- "Vaginularia_paradoxa" <- "4"
- "Humata repens" <- "12"

Add the following character scores:
Humata_repens, ribbon -> cordate, glands ? -> 1

(Date uncertain): Copied morph_qual_traits_2016-09-03_mod.csv from Moorea_ferns_morph_2016-09-03.xls and modified. Not sure when this happened, but I verified the changes with check_diff.R.

Changes to morph_qual_traits_2016-09-03_mod.csv from Moorea_ferns_morph_2016-09-03.xls as follows
(see check_diff.R for code):
- Delete the following rows (by OTU): "Austrobaileya_scandens", "Chloranthus_japonicus",
"Cycas_circinalis", "Diplazium_proliferum", "Elaphoglossum_florencei", "Ginkgo_biloba",
"Gnetum_gnemon", "Huperzia_phlegmaria", "Huperzia_ribourtii", "Huperzia_squarrosa", "Isoetes_melanopoda", "Lycopodiella_cernua", "Marsilea_polycarpa", "Ophioglossum_reticulatum",
"Pinus_radiata", "Prosaptia_subnuda", "Selaginella_apoda", "Selaginella_banksii", 
"Selaginella_laxa" (mostly outgroups)
- Change "key group" column by deleting commas
- Dicranopteris_linearis dissection from NA to 5
- Hymenophyllum_braithwaitei hairs from 0 to NA
- Hypolepis_sp1 hairs from ? to NA
- Lindsaea_propinqua hairs from ? to NA
- Sticherus_tahitensis hairs from ? to NA, dissection from NA to 5
