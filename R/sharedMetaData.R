#' Add Study to meta-data
#'
#' @param dt meta-data data.table
#' @export
#'
addStudy <- function(dt){
  dt$study_accession <- gsub("SUB\\d{6,7}\\.", "SDY", dt$participant_id)
  return(dt)
}

#' Add Arm Accession to meta-data
#'
#' @param dt meta-data data.table
#' @param gef gene expression files table from ImmuneSpace connection
#' @export
#'
addArmAccession <- function(dt, gef){
  dt$arm_accession <- gef$arm_accession[ match(dt$participant_id, gef$participant_id)]
  return(dt)
}

#' Add geBatchName for use in cross-study normalization
#'
#' @param dt meta-data data.table
#' @export
#'
addGeBatchName <- function(dt){
  # SDY1276 - Discovery/Male (ARM4428), Validation/Female (ARM4431)
  # SDY1264 - Trial1 (ARM4368), Trial2 (ARM4369)
  byCohort <- c("ARM4428", "ARM4431", "ARM4368", "ARM4369")

  # SDY180 - originally 3 groups with 2 vaccines per group, except Group 3
  # was saline and therefore not included
  byGroup <- c(ARM773 = "Group 1",
               ARM774 = "Group 1",
               ARM776 = "Group 2",
               ARM777 = "Group 2")

  dt$geBatchName <- apply(dt, 1, function(x){
    if(x[['arm_accession']] %in% byCohort){
      return(paste(x[['study_accession']], x[['arm_accession']], sep = "_"))
    }else if(x[['arm_accession']] %in% names(byGroup)){
      grpName <- byGroup[ match(x[['arm_accession']], names(byGroup))]
      return(paste(x[['study_accession']], grpName, sep = "_"))
    }else{
      return(x[['study_accession']])
    }
  })

  return(dt)
}

#' Add irpBatchName for use in responseCall generation
#'
#' @param dt meta-data data.table
#' @export
#'
addIrpBatchName <- function(dt){

  # Ensure that each vaccine within studies is looked at separately
  dt$irpBatchName <- paste(dt$study_accession, dt$vaccine, sep = "_")

  # For studies with significantly different cohorts over time
  byCohort <- c("ARM4428", "ARM4431", "ARM4368", "ARM4369")
  # SDY1276 - Discovery/Male (ARM4428), Validation/Female (ARM4431)
  # SDY1264 - Trial1 (ARM4368), Trial2 (ARM4369)

  dt$irpBatchName <- apply(dt, 1, function(x){
    if(x[['arm_accession']] %in% byCohort){
      return(paste(x[['irpBatchName']], x[['cohort']], sep = "_"))
    }else{
      return(x[['irpBatchName']])
    }
  })

  return(dt)
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
  dt$adjuvant <- vaccines$adjuvant[ armOrder ]
  return(dt)
}

#' Filter out cohorts with no vaccine data from custom mapping
#'
#' @param dt meta-data data.table
#' @export
#'
filterOutNoVaccineSamples <- function(dt){
  dt <- dt[ !is.na(dt$vaccine) & dt$vaccine_type != '', ]
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

#' Add fields from SDY1325 custom meta-data to shared meta-data
#'
#' @param dt meta-data data.table
#' @export
#'
addSDY1325metadata <- function(dt){
  sdy1325pids <- which(dt$participant_id %in% sdy1325_metadata$participant_id)
  matchOrder <- order(match(sdy1325_metadata$participant_id, dt$participant_id[sdy1325pids]))
  dt$age_reported[ sdy1325pids ] <- sdy1325_metadata$Age[ matchOrder ]
  dt$gender[ sdy1325pids ] <- sdy1325_metadata$Sex[ matchOrder ]
  return(dt)
}
