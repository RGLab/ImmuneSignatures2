outputDir = here::here("outputs", "2020_12_23_newAnno")
dataCacheDir = here::here("data_cache")
timestamp = "2020_12_23_newAnno_"

if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
if (!dir.exists(dataCacheDir)) dir.create(dataCacheDir, recursive = TRUE)

# rmarkdown::render(input = here::here("data-raw", "generate_base_eset.Rmd"),
#        output_file = here::here(outputDir, "generate_base_eset.html"),
#        params = list(
#          outputDir = outputDir,
#          dataCacheDir = dataCacheDir,
#          timestamp = timestamp
#        ))
# rmarkdown::render(input = here::here("data-raw", "create_final_esets.Rmd"),
#        output_file = here::here(outputDir, "create_final_esets.html"),
#        params = list(
#          outputDir = outputDir,
#          dataCacheDir = dataCacheDir,
#          timestamp = timestamp
#        ))
rmarkdown::render(input = here::here("vignettes", "pca_plots.Rmd"),
                  output_file = here::here(outputDir, "pca_plots.html"),
                  params = list(
                          outputDir = file.path(outputDir, "pca"),
                          dataCacheDir = dataCacheDir,
                          timestamp = timestamp
                  ))
