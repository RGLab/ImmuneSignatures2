#' Test gene expression matrices prior to pre-processing
#'
#' @param esets list of expressionSet objects
#' @export
#'
testExtractedGEData <- function(esets){

  chks <- list()
  chks$hasNames <- all(!is.null(names(esets)))

  naChk <- sapply(esets, function(x){
    em <- Biobase::exprs(x)
    return(any(is.na(em)))
  })

  chks$esetsWithNAvalues <- length(names(esets)[ naChk == TRUE ]) == 0

  dims <- lapply(esets, function(x){ dim(Biobase::exprs(x)) })

  # SDY1293 has fewest features - 19855 (subset from gsm files with 22207)
  dimChk <- unlist(sapply(dims, function(x){
    genes <- x[[1]] > 19500
    subs <- x[[2]] > 0
    return( genes & subs )
  }))

  chks$esetsWithBadDims <- length(names(esets)[ !dimChk ]) == 0

  return(chks)
}

#' Test gene expression meta-data prior to summarization by gene symbol
#'
#' @param geMetaData gene expression meta-data
#' @export
#'
testGEMetaDataPreSummarization <- function(geMetaData){
  results <- all(expectedGeMetaDataColumns %in% colnames(geMetaData))

  if(results){
    tmp <- geMetaData[ geMetaData$study_accession == "SDY1276", ]
    return(all(grepl("SDY1276", unique(tmp$matrix))))
  }else{
    return(results)
  }
}

#' Test gene expression matrix of all samples prior to cross-study normalization
#'
#' @param allGE gene expression matrix of all samples
#' @export
#'
testAllGEMatrixPreNorm <- function(allGE){
  chks <- list()

  chks$expectedNumberOfSubjects <- dim(allGE)[[2]] == 4869
  chks$expectedNumberOfGenes <- dim(allGE)[[1]] > 19500
  chks$completeCases <- sum(stats::complete.cases(allGE)) > 10000

  chks$importantGenesPresent <- all(c("MVP","ACTB") %in% allGE$rn)

  return(chks)
}

#' Test gene expression meta-data of all samples prior to cross-study normalization
#'
#' @param geMetaData gene expression meta-data of all samples
#' @export
#'
testAllGEMetaDataPreNorm <- function(geMetaData){
  chks <- list()

  tmp <- apply(geMetaData, 2, function(x){ all(is.na(x)) })
  tmp <- tmp[ !names(tmp) %in% c("exposure_material_reported", "exposure_process_preferred")]
  chks$naCols <- sum(tmp == TRUE) == 0

  expectedCellTypes <- c("Whole blood", "PBMC")
  chks$cellTypes <- length(setdiff(unique(geMetaData$cell_types), expectedCellTypes)) == 0

  yaleStudies <- c("SDY63", "SDY404", "SDY400", "SDY520", "SDY640")
  yaleTimepoints <- unique(geMetaData$time_post_last_vax[ geMetaData$study_accession %in% yaleStudies ])
  goodTimepoints <- c(0, 2, 4, 7, 28)
  chks$yaleTimepoints <- length(setdiff(yaleTimepoints, goodTimepoints)) == 0

  chks$sdy212biosampleRemoved <- !any(geMetaData$biosample_accession == "BS694717")

  chks$sdy1325day35Removed <- !any(geMetaData$time_post_last_vax[ geMetaData$study_accession == "SDY1325"] == 35)

  chks$noHourlyData <- !any(geMetaData$unit_post_last_vax == "Hours")

  chks$allMatricesPresent <- length(unique(geMetaData$matrix)) == 57

  chks$allAgesImputed <- all(geMetaData$age_imputed > 0 & !is.na(geMetaData$age_imputed))

  return(chks)
}

#' Test non-normalized base expression set for reasonable y-chromosome imputation and
#' removal of participants without baseline data
#'
#' @param eset non-normalized expression set with all participants
#' @export
#'
testNoNormEset <- function(eset){
  chks <- list()
  pd <- pData(eset)

  allPids <- unique(pd$participant_id)
  pidsWithBaseline <- unique(pd$participant_id[ pd$time_post_last_vax >= -7 & pd$time_post_last_vax <= 0 ])
  pidsToRm <- setdiff(allPids, pidsWithBaseline)
  chks$allPidsHaveBaseline <- length(pidsToRm) == 0

  chks$allPidsHaveImputedYchrom <- !any(is.na(pd$y_chrom_present) | pd$y_chrom_present == "NP")

  updatedGender <- pd[ grepl("male", pd$gender, ignore.case = TRUE), ]
  updatedGender$gender_imputed <- ifelse(updatedGender$y_chrom_present, "Male", "Female")
  updatedGender <- updatedGender[ updatedGender$gender != updatedGender$gender_imputed, ]
  pidsWithUpdatedGender <- length(unique(updatedGender$participant_id))
  chks$reasonableGenderImputation <- pidsWithUpdatedGender < 40

  problemSamples <- qualityControl.yChromPresentByMatrix(eset, returnObject = "probSamplesDT")
  chks$expectedYchromProblemSamples <- nrow(problemSamples) == 65

  chks$allStudiesPresent <- length(unique(pd$study_accession)) == 30

  return(chks)
}
