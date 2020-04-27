# Read in data for ELISA from internal source
library(data.table)
dt <- data.table::fread("data-raw/SDY1370_ELISA.csv", header = TRUE, skip = 1)
colsToKeep <- c("Case #", "Day post", "IgM OD-COV", "IgG OD-COV")
dt <- dt[, ..colsToKeep]

s2c <- data.table::fread("data-raw/SDY1370_subjectIdMap.csv")
s2c$arm_accession <- gsub("\\((ARM\\d{4}),.*", "\\1", s2c$Arm2Subject)

# Map original ids to ImmuneSpace ids
pidOrder <- match(dt$`Case #`, s2c$`User Defined ID`)
dt$participant_id <- paste0(s2c$Accession[ pidOrder ], ".1370")
dt[, c("participant_id",
       "arm_accession")
   :=
     list(paste0(s2c$Accession[ pidOrder ], ".1370"),
          s2c$arm_accession[ pidOrder])
   ]
dt[ , `Case #` := NULL ]

# update study_time_collected colname
setnames(dt, c("Day post", "IgG OD-COV", "IgM OD-COV"), c("study_time_collected", "IgG", "IgM"))

# dcast IgM and IgG to analyte col with value_reported as value
dt <- melt(dt,
           id.vars = c("participant_id", "study_time_collected"),
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
