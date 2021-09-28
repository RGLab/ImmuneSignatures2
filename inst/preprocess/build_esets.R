outputDir = here::here("outputs", "tmp")
dataCacheDir = here::here("data_cache", "2021_03_08")
timestamp = "2021_03_08_"

if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
if (!dir.exists(dataCacheDir)) dir.create(dataCacheDir, recursive = TRUE)

rmarkdown::render(input = here::here("data-raw", "generate_base_eset.Rmd"),
       output_file = here::here(outputDir, "generate_base_eset.html"),
       params = list(
         outputDir = outputDir,
         dataCacheDir = dataCacheDir,
         timestamp = timestamp
       ))
rmarkdown::render(input = here::here("data-raw", "create_final_esets.Rmd"),
       output_file = here::here(outputDir, "create_final_esets.html"),
       params = list(
         outputDir = outputDir,
         dataCacheDir = dataCacheDir,
         timestamp = timestamp
       ))
# rmarkdown::render(input = here::here("vignettes", "pca_plots.Rmd"),
#                   output_file = here::here(outputDir, "pca_plots.html"),
#                   params = list(
#                           outputDir = file.path(outputDir, "pca"),
#                           dataCacheDir = dataCacheDir,
#                           timestamp = timestamp
#                   ))

# rmarkdown::render(input = here::here("vignettes", "original_code", "IOF_RAPToR_age_imputation_all_functions_021121.Rmd"),
#                   output_file = here::here(outputDir, "RAPToR_all_functions.html"),
#                   output_format = "html_document",
#                   params = list(
#                     dataCacheDir = dataCacheDir,
#                     timestamp = timestamp
#                   ))
#
# rmarkdown::render(input = here::here("vignettes", "original_code", "IOF_RAPToR_age_imputation_best_functions_021121.Rmd"),
#                   output_file = here::here(outputDir, "RAPToR_best_functions.html"),
#                   output_format = "html_document",
#                   params = list(
#                     dataCacheDir = dataCacheDir,
#                     timestamp = timestamp
#                   ))
