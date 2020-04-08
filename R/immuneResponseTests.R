#' Ensure necessary columns are present in each applicable assay
#'
#' @param immdata nested assay data frame lists
#' @export
#'
testImmuneResponseData <- function(selectedImmdata){

      results <- lapply(selectedImmdata, function(dt){
        chks <- list()
        staticCols <- c('participant_id', 'study_accession', 'arm_accession',
                        'vaccine','pathogen','vaccine_type', 'adjuvant',
                        'ImmResp_baseline_value_MFC','ImmResp_baseline_timepoint_MFC',
                        'ImmResp_postVax_value_MFC','ImmResp_postVax_timepoint_MFC',
                        'maxStrain_RBA', 'maxStrain_MFC',
                        'ImmResp_baseline_value_RBA','ImmResp_baseline_timepoint_RBA',
                        'ImmResp_postVax_value_RBA','ImmResp_postVax_timepoint_RBA')

        chks$missingStaticCols <- length(setdiff(staticCols, colnames(dt))) == 0

        chks$missingDynamicCols <- length(grep("_p\\d+$", colnames(dt), value = TRUE)) >= 2

        chks$noMissingMFC <- all(!is.na(dt$MFC))

        chks$noBlankOrNAVaxType <- all(!is.na(dt$vaccine_type) & dt$vaccine_type != '')

        return(chks)
      })

  return(results)
}
