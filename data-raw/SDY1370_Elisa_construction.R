# Read in data for ELISA from internal source
library(data.table)
dt <- data.table::fread("data-raw/SDY1370_ELISA.csv", header = TRUE, skip = 1)
colsToKeep <- c("Case #", "Day post", "IgM OD-COV", "IgG OD-COV")
dt <- dt[, ..colsToKeep]

# Create map of author-given subject IDs to ImmuneSpace via GEO
library(ImmuneSpaceR)
library(GEOquery)
con <- CreateConnection("SDY1370", onTest = TRUE)
ge <- con$getDataset("gene_expression_files", original_view = TRUE)
ge <- ge[ ge$type == "PBMC" & ge$study_time_collected == 0, ]

sub2cohort <- apply(ge, 1, function(x){
  subject <- x[['participant_id']]
  gsm <- x[['geo_accession']]
  tmp <- getGEO(gsm)
  title <- tmp@header$title
  authorId <- gsub(" Day \\d{1,2}", "", gsub("PBMC ", "", title))
  return(c(subject, authorId))
})
s2c <- data.frame(t(sub2cohort))
colnames(s2c) <- c("subjectId", "authorId")

# Map original ids to ImmuneSpace ids
dt$participant_id <- s2c$subjectId[ match(dt$`Case #`, s2c$authorId)]
dt[ , `Case #` := NULL ]

# Missing people - 3020-6 and 3014-4
dt <- dt[ !is.na(participant_id)]

# update study_time_collected colname
setnames(dt, c("Day post", "IgG OD-COV", "IgM OD-COV"), c("study_time_collected", "IgG", "IgM"))

# dcast IgM and IgG to analyte col with value_reported as value
dt <- melt(dt,
           id.vars=c("participant_id", "study_time_collected"),
            measure.vars = c("IgG", "IgM"))
setnames(dt, c("variable", "value"), c("analyte", "value_reported"))

# Add necessary Columns
dt <- dt[, c("value_preferred",
             "unit_reported",
             "unit_preferred",
             "study_accession",
             "study_time_collected_unit",
             "experiment_accession",
             "comments")
         :=
           list(value_reported,
                "OD",
                "OD",
                "SDY1370",
                "Days",
                NA,
                NA)]


# Add arm_accession mapped from GE
dt$arm_accession <- ge$arm_accession[ match(dt$participant_id, ge$participant_id) ]

# Need to generate biosample and expsample accesions to be valid
tmp <- unique(paste(dt$participant_id, dt$study_time_collected))
rndAcc <- seq(999999 - length(tmp) + 1, 999999)
df <- data.frame(uid = tmp,
                 biosample_accession = paste0("BS", rndAcc),
                 expsample_accession = paste0("ES", rndAcc))
ord <- match(paste(dt$participant_id, dt$study_time_collected),
             df$uid)

dt[, c("biosample_accession",
       "expsample_accession")
   :=
     list( df$biosample_accession[ ord ],
           df$expsample_accession[ ord ])
   ]

# Set order
cn <- "participant_id, arm_accession, biosample_accession, expsample_accession, experiment_accession, study_accession, comments, analyte, study_time_collected, study_time_collected_unit, value_reported, value_preferred, unit_reported, unit_preferred"
cn <- strsplit(cn, ", ")[[1]]
setcolorder(dt, match(colnames(dt), cn))

# save out
sdy1370_elisa <- dt
usethis::use_data(sdy1370_elisa, overwrite = TRUE)
