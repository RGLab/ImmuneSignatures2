#' Ensure eset is as expected
#'
#' @param eset expressionSet object
#' @param expectResponse whether the eset has response call data
#' @param expectNormalization has the eset been normalized
#' @param ages integer vector of expected ages
#' @export
#'
testFinalEset <- function(eset, expectResponse, expectNormalization, ages){
  chks <- list(pdata = list(),
               exprs = list())

  # # expected number of subjects
  # expectedSamples <- list(withResponse = list(young = 2738,
  #                                             older = 1066),
  #                         noResponse = list(young = 3254,
  #                                           older = 1253))
  #
  # # Only for noNorm
  # expectedStudies <- list(withResponse = list(young = c("SDY1119, SDY1260, SDY1264, SDY1276, SDY1289, SDY1294, SDY1325, SDY1328, SDY1364, SDY1370, SDY1529, SDY180, SDY212, SDY224, SDY269, SDY270, SDY400, SDY404, SDY520, SDY56, SDY61, SDY63, SDY640, SDY80, SDY984"),
  #                                             older = c("SDY1119, SDY1328, SDY212, SDY400, SDY404, SDY520, SDY56, SDY63, SDY640, SDY67, SDY80, SDY984")
  #                                             ),
  #                         noResponse = list(young = c("SDY1119, SDY1260, SDY1264, SDY1276, SDY1289, SDY1291, SDY1293, SDY1294, SDY1325, SDY1328, SDY1364, SDY1370, SDY1373, SDY1529, SDY180, SDY212, SDY224, SDY269, SDY270, SDY400, SDY404, SDY520, SDY56, SDY61, SDY63, SDY640, SDY80, SDY984"),
  #                                           older = c("SDY1119, SDY1328, SDY1368, SDY212, SDY400, SDY404, SDY520, SDY56, SDY63, SDY640, SDY67, SDY80, SDY984")
  #                                           )
  #                         )

  noResponseCols <- c(expectedGeMetaDataColumns,
                      "gender_imputed",
                      "gender_imputed_timepoint",
                      "failedGenderQC")

  staticResponseCols <- c(
    'ImmResp_baseline_value_MFC','ImmResp_baseline_timepoint_MFC',
    'ImmResp_postVax_value_MFC','ImmResp_postVax_timepoint_MFC',
    'maxStrain_RBA', 'maxStrain_MFC',
    'ImmResp_baseline_value_RBA','ImmResp_baseline_timepoint_RBA',
    'ImmResp_postVax_value_RBA','ImmResp_postVax_timepoint_RBA'
  )

  expectedCols <- ifelse(expectResponse,
                         c(noResponseCols, staticResponseCols),
                         noResponseCols)

  expectedLevels <- paste0(c("low", "moderate", "high"), "Responder")

  # pData
  pd <- pData(eset)

  chks$pdata$allStaticColsPresent <- all(expectedCols %in% colnames(pd))
  chks$pdata$allUniqueUids <- sum(duplicated(pd$uid)) == 0

  uniquePids <- unique(pd$participant_id)
  pidsWithBaseline <- pd$participant_id[ pd$time_post_last_vax <= 0 | pd$time_post_last_vax >= -7 ]
  chks$pdata$allPidsHaveBaseline <- all(uniquePids %in% pidsWithBaseline)

  chks$pdata$expectedNumberOfSamples <- dim(pd)[[1]] > 1000

  chks$pdata$agesOk <- ifelse(length(ages) == 2,
                              all(pd$age_imputed >= ages[[1]] & pd$age_imputed < ages[[2]]),
                              all((pd$age_imputed >= ages[[1]] & pd$age_imputed < ages[[2]]) |
                                   pd$age_imputed >= ages[[3]] & pd$age_imputed < ages[[4]]))

  if(isTRUE(expectResponse)){
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

  # IS1 genes not expected in cross-study normalized esets
  if(!isTRUE(expectNormalization)){
    chks$exprs$IS1genes <- all(c("ACTB", "MVP") %in% unique(rownames(em)))
  }

  # Integration
  chks$namesMatch <- all.equal(colnames(em), pd$uid)

  return(chks)
}
