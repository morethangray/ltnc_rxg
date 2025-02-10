### A repository for the generalized mixed models used in: 

Valliere JM, Gray M, Ruygt J, Comendant T, Palladini M. (in review). Dry-season grazing enhances native plant diversity of oak savannas.

------------------------------------------------------------------------

### Analysis scripts

------------------------------------------------------------------------

### Data

#### .input/data

There are 6 data files within the data subfolder that are used in this analysis. Four .xlsx files contain the annual plant surveys conducted at Wantrup Preserve between 2019-2022.

#### wantrup_2019.xlsx

#### wantrup_2020.xlsx

#### wantrup_2021.xlsx

#### wantrup_2022.xlsx

Two files contain the abundance and richness estimates derived from the surveys.

#### abundance.csv

#### richness.csv

<br>

#### .input/lookup_tables

There are 8 data files within the lookup_files subfolder that are used in this analysis. They include:

#### model_terms.csv

| Column header | Data type | Description |
|------------------|---------------|----------------------------------------|
| abbr_term | Character | The abbreviated model term used in the R code. `treatment` indicates grazing treatment. `f_year` indicates survey year. `plot_type` indicates plot configuration (interior or exterior). `f_break` indicates the interval between survey periods (. `f_new` indicates a new plot (categorical yes, no). `f_one_yr` indicates a 1-year gap in grazing between surveys (categorical yes, no). `f_two_yr` indicates a 2-year gap in grazing between surveys (categorical yes, no). `grazer` indicates ruminant used (sheep or goat). |
| term | Character | Full name of the model term used to label results. |

<br>

#### .output/

------------------------------------------------------------------------

#### .R/

#### .R/functions
