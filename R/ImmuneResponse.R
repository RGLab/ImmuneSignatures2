#' Function generates a base dataset from ImmuneSpace
#'
#' @param assay assay name in ImmuneSpace
#' @param con ImmuneSpaceR connection object
#' @param studies list of ImmuneSpace studies to use in filtering
#' @export
#'
getImmuneResponseData <- function(assay, con, studies){
  dt <- con$getDataset(assayname, original_view = TRUE)
  dt <- dt[ dt$study_accession %in% studies, ]
}

#' Add fields from vaccine data to meta-data
#'
#' @param dt meta-data data.table
#' @param vaccineData vaccine data data.table
#' @export
#'
addVaccineFields <- function(dt, vaccineData){
  armOrder <- match(dt$arm_accession, vaccineData$arm_accession)
  dt$vaccine <- vaccineData$vaccine[ armOrder ]
  dt$vaccine_type <- vaccineData$vaccine_type[ armOrder ]
  dt$pathogen <- vaccineData$pathogen[ armOrder ]
  return(dt)
}

#' Filter out cohorts with no vaccine data from custom mapping
#'
#' @param dt meta-data data.table
#' @export
#'
filterOutNoVaccineSamples <- function(dt){
  dt <- dt[ !is.na(dt$vaccine), ]
}

#' Add fields from demographic data to meta-data
#'
#' @param dt meta-data data.table
#' @param demographicData demographic data data.table
#' @export
#'
addDemographicFields <- function(dt, demographicData){
  pidOrder <-  match(dt$participant_id, demographicData$participant_id)
  dt$age_reported <- demographicData$age_reported[ pidOrder ]
  dt$gender_reported <- demographicData$gender[ pidOrder ]
  return(dt)
}

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


#' Impute age for those with missing age values or ranges given by ImmPort
#'
#' @param dt meta-data data.table
#' @export
#'
imputeAge <- function(dt){

  fixAge <- function(study_accession, age_reported){
    if(study_accession %in% c("SDY1260", "SDY1264", "SDY1293")){
      # The following studies all have ranges of 18 to 45 years old described in their study meta-data.
      # Therefore, collaborators decided to use the mid-point of this range as the imputed age.
      return(30)
    }else if(study_accession == "SDY1370"){
      # The following study has a range of 18 to 55 years based on being a Phase 0/1 trial.
      # Therefore, collaborators decided to use the mid-point of this range as the imputed age.
      return(37)
    }else{
      return(age_reported)
    }
  }

  dt[, age_imputed := fixAge(study_accession, age_reported) ]

  impute.na <- function(x){ x[is.na(x)] <- mean(x, na.rm = TRUE); return(x) }

  for(study in unique(dt$study_accession)){
    allSmpls <- dt[ which(dt$study == study) ]
    noAgeSmpls <- allSmpls[ which(is.na(allSmpls$age_reported)) ]
    if ((length(noAgeSmpls) > 0) & (length(allSmpls) > length(noAgeSmpls))) {
      uidOrder <- match(allSmpls, dt$uid)
      dt$age_imputed[ uidOrder ] <- impute.na(dt$age_reported[uidOrder])
    }
  }

  return(dt)
}

