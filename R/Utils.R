#' Convert 'Hours'-based study times to 'Days'-based
#'
#' @param dt meta-data data.table
#' @import data.table
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

#' Get a backend table not easily accessible through ImmuneSpaceR
#'
#' @param con ImmuneSpaceR connection object
#' @param schemaName schema name in ImmuneSpace DB
#' @param queryName query name in ImmuneSpace DB
#' @param biosamples biosample_accession id vector
#' @import Rlabkey
#' @export
#'
getTable <- function(con, schemaName, queryName, biosamples = NULL){
  if(!is.null(biosamples)){
    bsFilter <- Rlabkey::makeFilter('biosample_accession',
                                    'IN',
                                    paste(biosamples, collapse = ";"))
    dt <- labkey.selectRows(baseUrl = con$config$labkey.url.base,
                            folderPath = con$config$labkey.url.path,
                            schemaName = schemaName,
                            queryName = queryName,
                            colNameOpt = "rname",
                            colFilterOpt = bsFilter)
  }else{
    dt <- labkey.selectRows(baseUrl = con$config$labkey.url.base,
                            folderPath = con$config$labkey.url.path,
                            schemaName = schemaName,
                            queryName = queryName,
                            colNameOpt = "rname")
  }
  return(dt)
}
