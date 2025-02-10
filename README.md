### A repository for the generalized mixed models used in:

Valliere JM, Gray M, Ruygt J, Comendant T, Palladini M. (in review). Dry-season grazing enhances native plant diversity of oak savannas.

------------------------------------------------------------------------

### Analysis scripts

**.R**

This repository has 4 R scripts used for the analysis. They will be in the R directory. These include:

**.R/1_setup.R**

**.R/2_model_selection_abundance.R**

**.R/3_model_selection_richness.R**

**.R/4_make_mm_and_contrast_tables.R**

<br>

**.R/functions**

This repository has 14 functions used for the analysis. They will be in the functions sub-directory. These include:

[Setup functions]{.underline}

-   **.R/functions/fxn_setup_dependencies.R**

-   **.R/functions/fxn_setup_file_paths.R**

-   **.R/functions/fxn_setup_lookup_tables.R**

[Data preparation functions]{.underline}

-   **.R/functions/fxn_prepare_input_data.R**

-   **.R/functions/fxn_calculate_metrics.R**

-   **.R/functions/fxn_load_rich_abun.R**

[Model selection functions]{.underline}

-   **.R/functions/fxn_identify_distribution.R**

-   **.R/functions/fxn_kable.R**

-   **.R/functions/fxn_lmer_fixed_effects.r**

-   **.R/functions/fxn_lmer_initial_model.r**

-   **.R/functions/fxn_lmer_model_review.r**

-   **.R/functions/fxn_lmer_random_effects.r**

-   **.R/functions/fxn_model_selection_rich.R**

-   fxn_summarize_models.R

------------------------------------------------------------------------

### **Data**

**.input/data**

There are 6 data files within the data subfolder that are used in this analysis.

Four .xlsx files contain the annual plant surveys conducted at Wantrup Preserve between 2019-2022. These were derived from the raw survey data and require pre-processing with `fxn_prepare_input_data`, which assigns consistent names to the column names, sheet names, and plant names (among other cleaning).

-   **wantrup_2019.xlsx**

-   **wantrup_2020.xlsx**

-   **wantrup_2021.xlsx**

-   **wantrup_2022.xlsx**

Two files contain the abundance and richness estimates derived from the surveys using `fxn_calculate_metrics`. Each file contains the estimates for each species subset (all native species, native forb species, and non-native species).

-   **abundance.csv**

-   **richness.csv:**

The script `fxn_load_rich_abun` is used to load the abundance and richness data into working memory as rich_abun.

<br>

**.input/lookup_tables**

There are 8 data files within the lookup_files subfolder that are used in this analysis. They include:

**model_terms.csv**

| Column header | Data type | Description |
|------------------|------------------|-------------------------------------|
| abbr_term | Character | The abbreviated model term used in the R code. `treatment` indicates grazing treatment. `f_year` indicates survey year. `plot_type` indicates plot configuration (interior or exterior). `f_break` indicates the interval between survey periods (. `f_new` indicates a new plot (categorical yes, no). `f_one_yr` indicates a 1-year gap in grazing between surveys (categorical yes, no). `f_two_yr` indicates a 2-year gap in grazing between surveys (categorical yes, no). `grazer` indicates ruminant used (sheep or goat). |
| term | Character | Full name of the model term used to label results. |

<br>

.output/

------------------------------------------------------------------------

**.R/**

**.R/functions**
