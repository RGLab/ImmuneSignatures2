#' Function generates a base dataset from ImmuneSpace
#'
#' @param assay assay name in ImmuneSpace
#' @param con ImmuneSpaceR connection object
#' @param studies list of ImmuneSpace studies to use in filtering
#' @export
#'
getImmuneResponseData <- function(assay, con, studies){
  dt <- con$getDataset(assay, original_view = TRUE)
  dt <- dt[ dt$study_accession %in% studies, ]
}

#' Add fields from vaccine data to meta-data
#'
#' @param dt meta-data data.table
#' @param vaccineData vaccine data data.table
#' @export
#'
addVaccineFields <- function(dt, vaccines){
  armOrder <- match(dt$arm_accession, vaccines$arm_accession)
  dt$vaccine <- vaccines$vaccine[ armOrder ]
  dt$vaccine_type <- vaccines$vaccine_type[ armOrder ]
  dt$pathogen <- vaccines$pathogen[ armOrder ]
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

  dt[, age_imputed := mapply(fixAge, study_accession, age_reported) ]

  fixNas <- function(ages){
    naLoc <- which(is.na(ages))
    if ((length(naLoc) > 0) & (length(ages) > length(naLoc))) {
      ages[naLoc] <- mean(ages, na.rm = TRUE)
    }
    return(ages)
  }

  dt[, age_imputed := fixNas(age_imputed), by = 'study_accession']

  return(dt)
}

#' Filter immdata list elements by age filter
#'
#' @param immdata list of assay data.table(s)
#' @param ageCutoff integer value for cutting age_imputed
#' @param isYoung boolean selecting young or older cohort
#' @export
#'
filterImmdataByAge <- function(immdata, ageCutoff, isYoung){
  if(isYoung){
    filteredImmdata <- lapply(immdata, function(dt){
      dt <- dt[dt$age_imputed < ageCutoff]
    })
  }else{
    filteredImmdata <- lapply(immdata, function(dt){
      dt <- dt[dt$age_imputed >= ageCutoff]
    })
  }
  return(filteredImmdata)
}
