# ImmuneSignatures2

This package provides the code, in the form of vignettes, to generate the expressionSet objects used in the ImmuneSignatures2 analysis by the Human Immunology Project Consortium as well as a vignette that reproduces the figures for the associated "Data Resource" manuscript.

### Vignettes:

_generate_base_eset.Rmd_
Constructs a single expressionSet object with transcriptomic data from over twenty public studies held in Immunespace (www.immunespace.org) and also collects immune response data for the same study participants into one list of dataframes.

_create_final_esets.Rmd_
Takes the single expressionSet object and immune response data created by `generate_base_eset.Rmd` and then generates 16 versions of this expressionSet based on participants' age, with and without response data and cross-study normalized or not. For more information on each of these procedures, please refer to the code and comments in the vignette.

- Age Group: young (18 to 50), old (60+), extended_old (50+), all
- Immune Response: with a response call (derived from HAI, NAb, or ELISA data) or without
- Cross-study normalized and batch-corrected: with and without

_data_resource_manuscript_figures.Rmd_
Reproduces all figures, including supplementary ones, for the associated manuscript using the expressionSet outputs from `create_final_esets.Rmd`.

_access_virtual_study_on_immunespace.Rmd_
The expressionSet objects created by the first two vignettes can be accessed using the ImmuneSpaceR API. This vignette demonstrates how to use the `CreateConnection("IS2")` method to connect to the virtual study and then access the assay data in various forms.

### Docker:

`/inst/docker/` contains a Dockerfile and associated files for building a docker container for reproducible analysis. This is used to create the docker image which all ImmuneSignatures2 code runs in on the ImmuneSpace servers. See the [README](/inst/docker/README.md) in that directory for more info.

### Notes:

- Original code by collaborators for the project is held in `vignettes/original_code` for reference
- The assay data used for the meta-analysis is effectively sandbboxed in ImmuneSpace as a "virtual" study and is not updated over time. When connecting to ImmuneSpace via the ImmuneSpaceR API using the "IS2" study id provides access to this snapshotted data and the annotation used at the time of manuscript publication.
- `vignettes/data_cache` is not tracked by git, but is used locally to hold the large outputs for the first two vignettes.
- `vignettes/outputs` is tracked by git and holds quality control plots from the first two vignettes as well as the rendered html for the `data_resource_manuscript_figures.Rmd`.
- There are some data used to generate the expressionSet objects (e.g. manually curated vaccine mapping for each study) that are not held in ImmuneSpace. All that data is found in the `data` directory and is created by vignettes found in the `data-raw` directory.
