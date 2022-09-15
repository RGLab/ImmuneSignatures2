vaccines <- data.table::fread("data-raw/vaccine_map.csv", header = TRUE)
vaccines <- vaccines[, V1 := NULL ]
usethis::use_data(vaccines, overwrite = TRUE)
