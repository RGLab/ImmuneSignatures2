#' Convert 'Hours'-based study times to 'Days'-based
#'
#' @param dt meta-data data.table
#' @export
#'
correctHrs <- function(dt){
  dt <- apply(dt, 1, function(row){
    if (row[["study_time_collected_unit"]] == "Hours"){
      row[["study_time_collected"]] <- as.numeric(row[["study_time_collected"]]) / 24
      row[["study_time_collected_unit"]] <- "Days"
    }
    return(row)
  })
  dt <- data.table(t(dt))
  dt$study_time_collected <- gsub(" ", "", dt$study_time_collected)
  dt$study_time_collected <- gsub("\\.00", "", dt$study_time_collected)
  return(dt)
}

#' Generate a unique id column
#'
#' @param dt meta-data data.table
#' @export
#'
createUniqueIdColumn <- function(dt){
  dt$uid <- paste(dt$participant_id,
                  dt$study_time_collected,
                  dt$study_time_collected_unit,
                  dt$biosample_accession,
                  sep = "_")
  return(dt)
}
