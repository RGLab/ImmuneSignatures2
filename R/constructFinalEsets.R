#' Add selected immune response data to final expressionSets
#'
#' @param eset expressionSets without response data
#' @param immdata selected immune response data
#' @export
#'
addResponseData <- function(eset, immdata){

    # Subsetting to participants with immune response data results
    # in there being genes without expression data. Must trim.
    eset <- eset[, eset$participant_id %in% immdata$participant_id ]
    eset <- removeAllNArows(eset)

    # TODO: Figure out how to handle cohorts!
    # rm cohort from immdata because it shows up as Young vs Old instead of arm name in eset
    immdata[, immdata$cohort := NULL]

    pd <- pData(eset)
    sharedCols <- intersect(colnames(pd), colnames(immdata))
    pdWithResponse <- merge(pd, immdata, by = sharedCols, all.x = TRUE)

    # Some participants drop due to insufficient data to make calls
    eset.withResponse <- eset[, eset$uid %in% pdWithResponse$uid]
    matchOrder <- order(match(pdWithResponse$uid, colnames(exprs(eset.withResponse))))
    pData(eset.withResponse) <- pdWithResponse[ matchOrder, ]

    return(eset.withResponse)
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
