vaccines <- data.table::fread("data-raw/vaccine_map_april_2020.csv", header = TRUE)
vaccines <- vaccines[, V1 := NULL ]
usethis::use_data(vaccines, overwrite = TRUE)
