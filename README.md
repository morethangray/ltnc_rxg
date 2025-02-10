### A repository for the generalized mixed models used in: 

Valliere JM, Gray M, Ruygt J, Comendant T, Palladini M. (in review). Dry-season grazing enhances native plant diversity of oak savannas.

------------------------------------------------------------------------

#### Analysis scripts

------------------------------------------------------------------------

#### Data

#### .input/data

#### .input/lookup_tables

##### model_terms.csv

| Column header | Data type | Description |
|------------------|---------------|----------------------------------------|
| abbr_term | Character | The abbreviated model term used in the R code. `treatment` indicates grazing treatment. `f_year` indicates survey year. `plot_type` indicates plot configuration (interior or exterior). `f_break` indicates the interval between survey periods (. `f_new` indicates a new plot (categorical yes, no). `f_one_yr` indicates a 1-year gap in grazing between surveys (categorical yes, no). `f_two_yr` indicates a 2-year gap in grazing between surveys (categorical yes, no). `grazer` indicates ruminant used (sheep or goat). |
| term | Character | Full name of the model term used to label results. |

#### .output/

------------------------------------------------------------------------

#### .R/

#### .R/functions
