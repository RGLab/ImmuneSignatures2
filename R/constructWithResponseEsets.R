addResponseData <- function(ageCohort, noResponseEsets, selectedImmdata){

    eset <- noResponseEsets[[ageCohort]]
    immdata <- selectedImmdata[[ageCohort]]

    eset <- eset[ , eset$participant_id %in% immdata$participant_id ]

    pd <- pData(eset)
    sharedCols <- intersect(colnames(pd), colnames(immdata))
    pData(eset) <- merge(pd, immdata, by = sharedCols)

    return(eset)
}
