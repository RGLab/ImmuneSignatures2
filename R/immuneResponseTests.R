#' Ensure necessary columns are present in each applicable assay
#'
#' @param immdata nested assay data frame lists
#' @export
#'
testImmuneResponseData <- function(immdata){
  results <- lapply(immdata, function(immdata_age_group){
    nms <- names(immdata_age_group)
    subResults <- lapply(nms, function(assay){
      if(assay == "elispot"){
        return(TRUE)
      }
      expectedCols <- c('vaccine','pathogen','vaccine_type',
                        'ImmResp_baseline_value_MFC','ImmResp_baseline_timepoint_MFC',
                        'ImmResp_postVax_value_MFC','ImmResp_postVax_timepoint_MFC',
                        'MFC_p30')
      if(assay %in% c('hai','neut_ab_titer')){
        expectedCols <- c(expectedCols,
                          'maxStrain_RBA', 'maxStrain_MFC', 'maxRBA_p30',
                          'ImmResp_baseline_value_RBA','ImmResp_baseline_timepoint_RBA',
                          'ImmResp_postVax_value_RBA','ImmResp_postVax_timepoint_RBA')
      }
      dt <- immdata_age_group[[assay]]
      res <- all(expectedCols %in% colnames(dt))
    })
    names(subResults) <- nms
    return(subResults)
  })

  return(results)
}
