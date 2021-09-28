# ImmuneSignatures2

This package provides the code used to generate the expressionSet objects used in the ImmuneSignatures2 analysis by the Human Immunology Project Consortium. 

### Preprocessing:

Preprocessing is run as a series of 3 R reports, implemented as an ETL in ImmuneSpace. These reports can be found in `inst/preprocess`: 

_pull_esets.Rmd_
Pulls 53 expressionsets from ImmuneSpace and saves them locally in an RDS object. 

_generate_base_eset.Rmd_
Constructs a single expressionSet object with transcriptomic data from over twenty public studies held in Immunespace (www.immunespace.org) and also collects immune response data for the same study participants into one list of dataframes.

_create_final_esets.Rmd_
Takes the single expressionSet object and immune response data created by `generate_base_eset.Rmd` and then generates 16 versions of this expressionSet based on participants' age, with and without response data and cross-study normalized or not. For more information on each of these procedures, please refer to the code and comments in the vignette.

- Age Group: young (18 to 50), old (60+), extended_old (50+), all
- Immune Response: with a response call (derived from HAI, NAb, or ELISA data) or without
- Cross-study normalized and batch-corrected: with and without

### Manuscript figures: 

Manuscript figures are reproduced on ImmuneSpace. The source for the Rmarkdown reports is saved in `inst/manuscript_figures`. These reports are run on ImmuneSpace in the docker container defined in this package, and outputs can be seen at ImmuneSpace.com/is2.url

### Supplementary preprocessing

Additional preprocessing code used in the Data Resource manuscript can be found in `inst/supplementary_preprocessing`. This includes code for additional figures not used in the manuscript, and some additional preprocessing steps.

### Package Data
Assay data is available on ImmuneSpace. Some additional datasets used in analysis and testing of data is saved as package data. Source code and data used to create these datasets can be found in `data-raw`. Package data includes:  
* `yChromGenes`: A vector of gene names from the Y-chromosome used in y-chromosome imputation. * `all_genesets`: A list of gene sets used in analysis. 
* `vaccines`: A data frame listing vaccines associated with each study present in the ImmuneSignatures2 data resource. 

### Docker:

`/inst/docker/` contains a Dockerfile and associated files for building a docker container for reproducible analysis. This is used to create the docker image which all ImmuneSignatures2 code runs in on the ImmuneSpace servers. See the [README](/inst/docker/README.md) in that directory for more info.

### Notes:

- Original code by collaborators for the project is held in `vignettes/original_code` for reference
- The assay data used for the meta-analysis is effectively sandbboxed in ImmuneSpace as a "virtual" study and is not updated over time. When connecting to ImmuneSpace via the ImmuneSpaceR API using the "IS2" study id provides access to this snapshotted data and the annotation used at the time of manuscript publication.
- `vignettes/data_cache` is not tracked by git, but is used locally to hold the large outputs for the first two vignettes.
- `vignettes/outputs` is tracked by git and holds quality control plots from the first two vignettes as well as the rendered html for the `data_resource_manuscript_figures.Rmd`.
- There are some data used to generate the expressionSet objects (e.g. manually curated vaccine mapping for each study) that are not held in ImmuneSpace. All that data is found in the `data` directory and is created by vignettes found in the `data-raw` directory.
