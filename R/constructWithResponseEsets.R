#' Add selected immune response data to final expressionSets
#'
#' @param ageCohort name of age cohort
#' @param noResponseEsets expressionSets without response data
#' @param selectedImmdata selected immune response data
#' @export
#'
addResponseData <- function(ageCohort, noResponseEsets, selectedImmdata){

    eset <- noResponseEsets[[ageCohort]]
    immdata <- selectedImmdata[[ageCohort]]

    eset <- eset[ , eset$participant_id %in% immdata$participant_id ]
    eset <- removeAllNArows(eset)

    pd <- pData(eset)
    sharedCols <- intersect(colnames(pd), colnames(immdata))
    pd <- merge(pd, immdata, by = sharedCols)

    pData(eset) <- pd[ order(match(pd$uid, colnames(exprs(eset)))), ]

    return(eset)
}
