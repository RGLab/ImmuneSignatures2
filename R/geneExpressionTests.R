#' Test gene expression matrices prior to pre-processing
#'
#' @param esets list of expressionSet objects
#' @import Biobase
#' @export
#'
testExtractedGEData <- function(esets){

  naChk <- sapply(esets, function(x){
    em <- Biobase::exprs(x)
    return(any(is.na(em)))
  })

  esetsWithNAvalues <- names(esets)[ naChk == TRUE ]

  dims <- lapply(eset, function(x){ dim(Biobase::exprs(x)) })

  # SDY1293 has fewest features - 19855 (subset from gsm files with 22207)
  dimChk <- lapply(dims, function(x){
    genes <- x[[1]] > 19500
    subs <- x[[2]] > 0
    return( genes & subs )
  })

  esetsWithBadDims <- names(esets)[ !dimChk]

  results <- list(NAvalues = esetsWithNAvalues,
                  badDims = esetsWithBadDims)

  return(results)
}
