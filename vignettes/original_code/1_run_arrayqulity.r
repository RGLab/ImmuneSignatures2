# Samples flagged as outliers by 3 / 3 quality scores will be removed from downstream analysis in the virtual study.
# This script outputs a csv file with the results with a header containing each of the QC metrics, and each array as rows.
# If the array i passed QC test j the qc test, the output is NA. If array i failed QC test j, the output has a 1.
# For each cohort n=9 summary plots and a detailed HTML report are created
# QC plot containing directory is named by the SDY accession number.
suppressMessages(library(tidyverse))
suppressMessages(library(arrayQualityMetrics))
suppressMessages(library(here))
suppressMessages(library(rlist))
options(stringsAsFactors = FALSE)
here()

# setup datapath
datapath = here("QC/generated_data/")
figpath = here("QC/figures/")
dir.create(datapath, recursive = TRUE); dir.create(figpath, recursive = TRUE)

# load virtual study
es = readRDS(file = here("data/2020_08_10_all_noNorm_withResponse_eset.rds"))

# index over unique SDY
stopifnot(!is_null(es$study_accession))
sdy = es$study_accession %>% unique

# build QC table
qc.list =  list()
for (i in 1:length(sdy)) {

  # make directory for output for sdy i
  stopifnot(!is.null(datapath))
  sdy_path = paste0(datapath, sdy[i]) ; dir.create(sdy_path)
  print(paste0(i, " of ", length(sdy),  " writing qc to ",  sdy_path))

  # subset to SDY acession
  array = es[ , es$study_accession == sdy[i]]

  # run AQM on sdy and save outlier table.
  qc = arrayQualityMetrics(expressionset = array,
                           outdir = sdy_path,
                           intgroup = c("matrix", "time_post_last_vax"),
                           force = TRUE,
                           spatial = FALSE,
                           do.logtransform = FALSE)
  if(!is.null(qc$arrayTable)){
    qc.table = qc$arrayTable
    qc.list[[i]] = qc.table
    } else {
      warning(paste0("no report table for ", sdy[i] ))
    }
}
# save list object
saveRDS(qc.list, file =paste0(datapath, "qclist.rds"))
