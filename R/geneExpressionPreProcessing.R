#' Add fields from immune exposure data to meta-data
#'
#' @param dt meta-data data.table
#' @param immuneExposureData immune exposure data.table
#' @export
#'
addImmuneExposureFields <- function(dt, immuneExposureData){
  pidOrder <-  match(dt$participant_id, expo$participant_id)
  dt$exposure_material <- immuneExposureData$`Exposure Material Reported`[ pidOrder ]
  dt$exposure_process <- immuneExposureData$`Exposure Process Preferred`[ pidOrder ]
  return(dt)
}
