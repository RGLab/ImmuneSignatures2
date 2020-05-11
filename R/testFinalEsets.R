#' Ensure eset is as expected
#'
#' @param eset expressionSet object
#' @param hasResponse whether the eset has response call data
#' @param ageCohort  whether the eset is for younger cohort
#' @param ageCutoffs the age cutoff point for the cohorts
#' @export
#'
testFinalEset <- function(eset, hasResponse, ageCohort, ageCutoffs){
  chks <- list(pdata = list(),
               exprs = list())

  # expected number of subjects
  expectedSamples <- list(withResponse = list(young = 2533,
                                              older = 931),
                          noResponse = list(young = 3045,
                                            older = 1111))

  expectedStudies <- list(withResponse = list(young = c("SDY269, SDY61, SDY270, SDY63, SDY224, SDY404, SDY400, SDY212, SDY56, SDY520, SDY640, SDY1119, SDY80, SDY180, SDY1289, SDY1276, SDY1294, SDY1264, SDY1325, SDY984, SDY1260, SDY67"),
                                              older = c("SDY63, SDY404, SDY400, SDY212, SDY56, SDY520, SDY640, SDY1119, SDY80, SDY984, SDY67")
                                              ),
                          noResponse = list(young = c("SDY520, SDY640, SDY80, SDY180, SDY1294, SDY1119, SDY1289, SDY1370, SDY1368, SDY67, SDY224, SDY212, SDY270, SDY1373, SDY1364, SDY1325, SDY1291, SDY1293, SDY1276, SDY1264, SDY1260, SDY984, SDY61, SDY56, SDY63, SDY404, SDY400, SDY269"),
                                            older = c("SDY520, SDY640, SDY80, SDY1119, SDY1368, SDY67, SDY212, SDY984, SDY56, SDY63, SDY404, SDY400")
                                            )
                          )

  noResponseCols <- c(expectedGeMetaDataColumns, "gender_imputed")

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

  chks$agesOk <- ifelse(ageCohort == "young",
                        all( pd$age_imputed < ageCutoffs[[1]] ),
                        all( pd$age_imputed >= ageCutoffs[[2]] )
                        )

  if(responseStatus == "withResponse"){
    chks$pdata$dynamicColsPresent <- any(grepl("^(MFC|maxRBA)_p\\d", colnames(pd)))

    mfcDiscretizedCols <- grep("^MFC_p\\d", colnames(pd))
    chks$pdata$noNAsInDiscretizedMFC <- !any(is.na(pd[[mfcDiscretizedCols[[1]] ]]))

    chks$pdata$expectedMFCLevels <- all(levels(pd[[mfcDiscretizedCols[[1]]]]) == expectedLevels)
  }

  # exprs
  em <- Biobase::exprs(eset)

  chks$exprs$genesWithCompleteCases <- sum(complete.cases(em)) > 10000
  chks$exprs$noMissingGeneNames <- all(!is.na(rownames(em)) & rownames(em) != "")

  incompleteRows <- apply(em, 1, function(x){ all(is.na(x)) })
  chks$exprs$noIncompleteRows <- sum(incompleteRows) == 0

  chks$IS1genes <- all(c("ACTB", "MVP") %in% unique(rownames(em)))

  # Integration
  chks$namesMatch <- all.equal(colnames(em), pd$uid)

  return(chks)
}
