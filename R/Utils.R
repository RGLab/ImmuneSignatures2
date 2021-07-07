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
#' @param ... additional arguments passed to labkey.selectRows.
#' @import Rlabkey
#' @export
#'
getTable <- function(con, schemaName, queryName, ...){
  dt <- labkey.selectRows(baseUrl = con$config$labkey.url.base,
                          folderPath = con$config$labkey.url.path,
                          schemaName = schemaName,
                          queryName = queryName,
                          colNameOpt = "rname",
                          ...)
  return(dt)
}


#' Write a log of processing date to a csv
#'
#' @param metadata_path The metadata csv file
#' @param task_name the name of the task
#' @export
write_processing_metadata <- function(metadata_path,
                                      task_name) {

  if ( file.exists(metadata_path) ) {
    processing_date <- fread(metadata_path)
  } else {
    processing_date <- data.table(task = task_name)
  }

  if ( task_name %in% processing_date$task ) {
    processing_date[task == task_name,
                    `:=`(
                      date = strftime(Sys.time(), "%Y_%m_%d", tz = "US/Pacific"),
                      ImmuneSignatures2_version = as.character(packageVersion("ImmuneSignatures2"))
                    )]
  } else {
    processing_date <- rbind(processing_date,
                             data.table(
                               task = task_name,
                               date = strftime(Sys.time(), "%Y_%m_%d", tz = "US/Pacific"),
                               ImmuneSignatures2_version = as.character(packageVersion("ImmuneSignatures2"))
                             ))
  }

  fwrite(processing_date, metadata_path)
}

#' Write metadata about data to a csv
#'
#' @param metadata_path The metadata csv file
#' @param dataset_name Name of the dataset
#' @param data_path Path to the dataset
#' @param include_counts Include count of subjects, samples and features?
#' if \code{TRUE}, \code{data} must not be \code{NULL}.
#' @param data the data
#' @export
write_data_metadata <- function(metadata_path,
                                dataset_name,
                                data_path,
                                data = NULL,
                                include_counts = FALSE) {


  if ( file.exists(metadata_path) ) {
    metadata <- fread(metadata_path)
  } else {
    metadata <- data.table(dataset = dataset_name)
  }

  if ( dataset_name %in% metadata$dataset ) {
    metadata[dataset == dataset_name,
             `:=`(
               path = data_path,
               date = strftime(Sys.time(), "%Y_%m_%d", tz = "US/Pacific"),
               ImmuneSignatures2_version = as.character(packageVersion("ImmuneSignatures2"))
             )]
  } else {
    metadata <- rbind(metadata,
                      data.table(
                        dataset = dataset_name,
                        path = data_path,
                        date = strftime(Sys.time(), "%Y_%m_%d", tz = "US/Pacific"),
                        ImmuneSignatures2_version = as.character(packageVersion("ImmuneSignatures2"))
                      ),
                      fill = TRUE)
  }


  if (include_counts & !is.null(data)) {
    if ( class(data) == "ExpressionSet" ) {
      metadata[dataset == dataset_name,
               `:=`(
                 subjects = length(unique(data$participant_id)),
                 samples = dim(data)["Samples"],
                 features = dim(data)["Features"]
               )]
    } else {
      metadata[dataset == dataset_name,
               `:=`(
                 subjects = length(unique(data$participant_id)),
                 samples = dim(data)[1],
                 featrues = dim(data)[2]
               )]
    }
  }
  fwrite(metadata, metadata_path)
}
