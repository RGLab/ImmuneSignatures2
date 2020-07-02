# Using Fig1 Data from Paper, creating mapping to biosample for merging into sharedMetaData
library(data.table)
dt <- fread("data-raw/SDY1325_Fig1_Data.csv")

library(ImmuneSpaceR)
con <- CreateConnection("SDY1325", onTest = TRUE)
gef <- con$getDataset("gene_expression_files", original_view = TRUE)
gef <- gef[ gef$study_time_collected == 0, ]
gef$authorId <- sapply(gef$geo_accession, function(gsm){
  rec <- GEOquery::getGEO(gsm)
  title <- rec@header$title
  authorId <- gsub("Adult (\\d{1,2}).+", "\\1", title)
})
gef$authorId <- as.integer(gef$authorId)

metadata <- merge(gef, dt, by.x = "authorId", by.y = "Participant")
colsToKeep <- c("participant_id", "Age", "Sex")
sdy1325_metadata <- metadata[ , ..colsToKeep ]
usethis::use_data(sdy1325_metadata, overwrite = TRUE)
