# ImmuneSignatures2

[![R-CMD-check](https://github.com/RGLab/ImmuneSignatures2/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/ImmuneSignatures2/actions)
[![DOI](https://zenodo.org/badge/252603828.svg)](https://zenodo.org/badge/latestdoi/252603828)
[![docker](https://github.com/RGLab/ImmuneSignatures2/actions/workflows/docker-build.yaml/badge.svg)](https://hub.docker.com/r/rglab/immunesignatures2)

This package provides the code used to generate the expressionSet objects used in the ImmuneSignatures2 analysis by the Human Immunology Project Consortium. Data, as well as replication of all figures, is available on [ImmuneSpace](www.immunespace.org/is2.url)

### Preprocessing:

Preprocessing is run as a series of 3 R reports, implemented as an ETL in ImmuneSpace. These reports can be found in `inst/preprocess`: 

[_pull_esets.Rmd_](/inst/preprocess/pull_esets.Rmd) 
Pulls 53 expressionsets from ImmuneSpace and saves them locally in an RDS object. 
These expressionsets have already been run through the basic ImmuneSpace preprocessing pipeline, and run through some basic normalization steps. The code used for the ImmuneSpace preprocessing is available [here](https://github.com/RGLab/LabKeyModules/blob/b21d1eaf1e67a209ec376ab6ee1a06e04c35aa41/HIPCMatrix/pipeline/tasks/runCreateMx.R#L576), where "NA" plaform corresponds to RNA-seq data. 

[_generate_base_eset.Rmd_](/inst/preprocess/generate_base_eset.Rmd) 
Constructs a single expressionSet object with transcriptomic data from over twenty public studies held in Immunespace (www.immunespace.org) and also collects immune response data for the same study participants into one list of dataframes.

[_create_final_esets.Rmd_](/inst/preprocess/create_final_esets.Rmd) 
Takes the single expressionSet object and immune response data created by `generate_base_eset.Rmd` and then generates 16 versions of this expressionSet based on participants' age, with and without response data and cross-study normalized or not. For more information on each of these procedures, please refer to the code and comments in the vignette.

- Age Group: young (18 to 50), old (60+), extended_old (50+), all
- Immune Response: with a response call (derived from HAI, NAb, or ELISA data) or without
- Cross-study normalized and batch-corrected: with and without

### Manuscript figures: 

Manuscript figures are reproduced on ImmuneSpace. The source for the Rmarkdown reports is saved in `inst/manuscript_figures`. These reports are run on ImmuneSpace in the docker container defined in this package, and outputs can be seen on the [ImmuneSignatures page](https://www.immunespace.org/is2.url) on ImmuneSpace. 

### Supplementary preprocessing

Additional preprocessing code used in the Data Resource manuscript can be found in `inst/supplementary_preprocessing`. This includes code for additional figures not used in the manuscript, and some additional preprocessing steps.

### Package Data
Assay data is available on ImmuneSpace. Some additional datasets used in analysis and testing of data is saved as package data. Source code and data used to create these datasets can be found in `data-raw`. Package data includes:   
* `yChromGenes`: A vector of gene names from the Y-chromosome used in y-chromosome imputation.  
* `all_genesets`: A list of gene sets used in analysis.  
* `vaccines`: A data frame listing vaccines associated with each study present in the ImmuneSignatures2 data resource.  

### Docker:

`/inst/docker/` contains a Dockerfile and associated files for building a docker container for reproducible analysis. This is used to create the docker image which all ImmuneSignatures2 code runs in on the ImmuneSpace servers. The docker image is available on [Dockerhub](https://hub.docker.com/r/rglab/immunesignatures2). See the [README](/inst/docker/README.md) in that directory for more info.

### Data Access

The assay data used for the meta-analysis is effectively sandbboxed in ImmuneSpace as a "virtual" study and is not updated over time. When connecting to ImmuneSpace via the ImmuneSpaceR API using the "IS2" study id provides access to this snapshotted data and the annotation used at the time of manuscript publication.
