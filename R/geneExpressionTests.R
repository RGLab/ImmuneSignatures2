#' Test gene expression matrices prior to pre-processing
#'
#' @param esets list of expressionSet objects
#' @import Biobase
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

  expectedCols <- c("uid",
                    "participant_id",
                    "biosample_accession",
                    "study_time_collected",
                    "study_time_collected_unit",
                    "time_post_last_vax",
                    "unit_post_last_vax",
                    "age_reported",
                    "age_imputed",
                    "gender",
                    "race",
                    "ethnicity",
                    "exposure_material",
                    "exposure_process",
                    "matrix",
                    "gsm",
                    "study",
                    "study2",
                    "Hispanic",
                    "White",
                    "Asian",
                    "Black",
                    "cell_type",
                    "cohort",
                    "fas_name",
                    "vendor"
  )

  results <- all(expectedCols %in% colnames(geMetaData))

  return(results)
}

#' Test gene expression matrix of all samples prior to cross-study normalization
#'
#' @param allGE gene expression matrix of all samples
#' @export
#'
testAllGEMatrixPreNorm <- function(allGE){
  chks <- list()

  chks$expectedNumberOfSubjects <- dim(allGE)[[2]] == 4354
  chks$expectedNumberOfGenes <- dim(allGE)[[1]] > 19500
  chks$completeCases <- sum(complete.cases(allGE)) > 10000

  return(chks)
}

#' Test gene expression meta-data of all samples prior to cross-study normalization
#'
#' @param geMetaData gene expression meta-data of all samples
#' @export
#'
testAllGEMetaDataPreNorm <- function(geMetaData){
  chks <- list()

  tmp <- apply(geMetaData, 2, function(x){ all(is.na(x))})
  tmp <- tmp[ !names(tmp) %in% c("exposure_material_reported", "exposure_process_preferred")]
  chks$naCols <- sum(tmp == TRUE) == 0

  sdy180cellTypes <- unique(geMetaData$cell_type[ geMetaData$study_accession == "SDY180"])
  sdy80cellTypes <- unique(geMetaData$cell_type[ geMetaData$study_accession == "SDY80"])
  chks$cellTypes <- sdy180cellTypes == "Whole blood" & sdy80cellTypes == "PBMC"

  yaleStudies <- c("SDY63", "SDY404", "SDY400", "SDY520", "SDY640")
  yaleTimepoints <- unique(geMetaData$time_post_last_vax[ geMetaData$study_accession %in% yaleStudies ])
  badTimepoints <- c(3, 5, 8, 9, 24)
  chks$yaleTimepoints <- !any(badTimepoints %in% as.numeric(yaleTimepoints))

  chks$sdy212biosampleRemoved <- !any(geMetaData$biosample_accession == "BS694717")

  chks$sdy1325day35Removed <- !any(geMetaData$time_post_last_vax[ geMetaData$study_accession == "SDY1325"] == 35)

  chks$noHourlyData <- !any(geMetaData$unit_post_last_vax == "Hours")

  chks$allMatricesPresent <- length(unique(geMetaData$matrix)) == 55

  chks$allAgesImputed <- all(geMetaData$age_imputed > 0 & !is.na(geMetaData$age_imputed))

  return(chks)
}

