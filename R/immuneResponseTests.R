#' Ensure necessary columns are present in each applicable assay
#'
#' @param immdata nested assay data frame lists
#' @export
#'
testImmuneResponseData <- function(selectedImmdata){

      expectedStudies <- list(young = c("SDY269, SDY61, SDY270, SDY63, SDY224, SDY404, SDY400, SDY212, SDY56, SDY520, SDY640, SDY1119, SDY80, SDY180, SDY1289, SDY1276, SDY1294, SDY1264, SDY1325, SDY984, SDY1260, SDY1364, SDY1328, SDY1370"),
                              older = c("SDY63, SDY404, SDY400, SDY212, SDY56, SDY520, SDY640, SDY1119, SDY80, SDY984, SDY67"))

      results <- lapply(names(selectedImmdata), function(ageCohort){

        chks <- list()

        dt <- selectedImmdata[[ageCohort]]
        studies <- strsplit(expectedStudies[[ageCohort]], ", ")[[1]]

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

        chks$hasExpectedStudies <- all(studies %in% unique(dt$study_accession))

        return(chks)
      })

  return(results)
}
