#' Ensure eset is as expected
#'
#' @param eset expressionSet object
#' @param hasResponse whether the eset has response call data
#' @param isYoung whether the eset is for younger cohort
#' @param ageCutff the age cutoff point for the cohorts
#' @export
#'
testFinalEset <- function(eset, responseStatus, ageCohort, ageCutoff){
  chks <- list(pdata = list(),
               exprs = list())

  # expected number of subjects
  expectedSamples <- list(withResponse = list(young = 2474,
                                           older = 787),
                       noResponse = list(young = 3380,
                                         older = 973))

  expectedGenes <- list(older = 15018,
                        young = 10056)

  noResponseCols <- c(
    'participant_id', 'study_accession', 'arm_accession', 'uid', 'cohort',
    'biosample_accession', 'time_post_last_vax', 'unit_post_last_vax',
    'age_imputed', 'race', 'gender_imputer', 'ethnicity',
    'matrix', 'gsm', 'cell_type', 'featureSetName', 'featureSetVendor',
    'vaccine','pathogen','vaccine_type','adjuvant'
  )

  staticResponseCols <- c(
    'ImmResp_baseline_value_MFC','ImmResp_baseline_timepoint_MFC',
    'ImmResp_postVax_value_MFC','ImmResp_postVax_timepoint_MFC',
    'maxStrain_RBA', 'maxStrain_MFC',
    'ImmResp_baseline_value_RBA','ImmResp_baseline_timepoint_RBA',
    'ImmResp_postVax_value_RBA','ImmResp_postVax_timepoint_RBA'
  )

  expectedCols <- ifelse(responseStatus == "noResponse",
                         noResponseCols,
                         c(noResponseCols, staticResponseCols))

  expectedLevels <- paste0(c("low", "moderate", "high"), "Responder")

  # pData
  pd <- pData(eset)

  chks$pdata$allStaticColsPresent <- all(expectedCols %in% colnames(pd))
  chks$pdata$allUniqueUids <- sum(duplicated(pd$uid)) == 0

  uniquePids <- unique(pd$participant_id)
  pidsWithBaseline <- pd$participant_id[ pd$time_post_last_vax <= 0 | pd$time_post_last_vax >= -7 ]
  chks$pdata$allPidsHaveBaseline <- all(uniquePids %in% pidsWithBaseline)

  chks$pdata$expectedNumberOfSamples <- dim(pd)[[1]] == expectedSamples[[responseStatus]][[ageCohort]]

  if(responseStatus == "withResponse"){
    chks$pdata$dynamicColsPresent <- any(grepl("^(MFC|maxRBA)_p\\d", colnames(pd)))

    mfcDiscretizedCols <- grep("^MFC_p\\d", colnames(pd))
    chks$pdata$noNAsInDiscretizedMFC <- !any(is.na(pd[[mfcDiscretizedCols[[1]] ]]))

    chks$pdata$expectedMFCLevels <- all(levels(pd[[mfcDiscretizedCols[[1]]]]) == expectedLevels)
  }

  # exprs
  em <- Biobase::exprs(eset)

  chks$exprs$genesWithCompleteCases <- sum(complete.cases(em)) == expectedGenes[[ageCohort]]
  chks$exprs$noMissingGeneNames <- all(!is.na(rownames(em)) & rownames(em) != "")

  incompleteRows <- apply(em, 1, function(x){ all(is.na(x)) })
  chks$exprs$noIncompleteRows <- sum(incompleteRows) == 0

  return(chks)
}
