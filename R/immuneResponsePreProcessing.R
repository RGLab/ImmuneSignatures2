#' Function generates a base dataset from ImmuneSpace
#'
#' @param assay assay name in ImmuneSpace
#' @param con ImmuneSpaceR connection object
#' @param studies list of ImmuneSpace studies to use in filtering
#' @export
#'
getImmuneResponseData <- function(assay, con, studies){
  dt <- con$getDataset(assay, original_view = TRUE)
  dt <- dt[ dt$study_accession %in% studies, ]
}

#' Filter immdata list elements by age filter
#'
#' @param immdata list of assay data.table(s)
#' @param ageCutoff integer value for cutting age_imputed
#' @param isYoung boolean selecting young or older cohort
#' @export
#'
filterImmdataByAge <- function(immdata, ageCutoffs, isYoung){
  if(isYoung){
    filteredImmdata <- lapply(immdata, function(dt){
      dt <- dt[dt$age_imputed < ageCutoffs[[1]]]
    })
  }else{
    filteredImmdata <- lapply(immdata, function(dt){
      dt <- dt[dt$age_imputed >= ageCutoffs[[2]]]
    })
  }
  return(filteredImmdata)
}

#' Generate a single response call data.table from multiple assays
#'
#' @param immdata_age_group list of assay data.table(s)
#' @export
#'
selectResponsesToUse <- function(immdata_age_group){

  tmp <- rbindlist(immdata_age_group, fill = TRUE)

  colsToUse <- c(
    "participant_id", "study_accession", "arm_accession", "cohort",
    "race", "ethnicity", "gender", "age_imputed",
    "vaccine", "vaccine_type", "pathogen", "assay", "adjuvant",
    "MFC", "maxRBA",
    "maxStrain_MFC", "maxStrain_RBA",
    "ImmResp_baseline_value_MFC", "ImmResp_baseline_timepoint_MFC",
    "ImmResp_postVax_value_MFC",  "ImmResp_postVax_timepoint_MFC",
    "ImmResp_baseline_value_RBA", "ImmResp_baseline_timepoint_RBA",
    "ImmResp_postVax_value_RBA",  "ImmResp_postVax_timepoint_RBA"
  )

  dynamicCols <- grep("_p\\d+$", colnames(tmp), value = TRUE)
  colsToUse <- c(colsToUse, dynamicCols)

  tmp <- tmp[, ..colsToUse]
  tmp <- tmp[ !duplicated(tmp) ] # after rm biosample id and other cols, mostly dupes
  tmp <- tmp[ , numAssays := .N, by = participant_id ]

  # select assay to use
  singleAssay <- tmp[ numAssays == 1 ]
  needFiltering <- tmp[ numAssays > 1 ]
  haiSubs <- needFiltering[ assay == "hai" ]
  nabSubs <- needFiltering[ assay == "neut_ab_titer" &
                            !participant_id %in% unique(haiSubs$participant_id) ]

  res <- rbindlist(list(singleAssay, haiSubs, nabSubs))

  if( length(res$participant_id) != length(unique(res$participant_id)) ){
    stop("filtering not done correctly")
  }

  return(res)
}
