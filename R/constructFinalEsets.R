#' Add selected immune response data to final expressionSets
#'
#' @param eset expressionSets without response data
#' @param immdata selected immune response data
#' @export
#'
addResponseData <- function(eset, immdata){

    eset <- eset[ , eset$participant_id %in% immdata$participant_id ]
    eset <- removeAllNArows(eset)

    pd <- pData(eset)
    sharedCols <- intersect(colnames(pd), colnames(immdata))
    pd <- merge(pd, immdata, by = sharedCols)

    pData(eset) <- pd[ order(match(pd$uid, colnames(exprs(eset)))), ]

    return(eset)
}

#' Filter expressionSet by age cutoffs
#'
#' @param noNormEset expressionSets without response data
#' @param ages allowed ages, either one or two sets of low and high points
#' @export
#'
filterEsetByAge <- function(noNormEset, ages){
    if(length(ages) == 2){
        eset <- noNormEset[ , noNormEset$age_imputed >= ages[[1]] &
                                noNormEset$age_imputed < ages[[2]] ]
    }else{
        eset <- noNormEset[ , (noNormEset$age_imputed >= ages[[1]] &
                                   noNormEset$age_imputed < ages[[2]]) |
                                (noNormEset$age_imputed >= ages[[3]] &
                                     noNormEset$age_imputed < ages[[4]]) ]
    }
    return(eset)
}
