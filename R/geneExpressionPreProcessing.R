#' #' Add fields from immune exposure data to meta-data
#' #'
#' #' @param dt meta-data data.table
#' #' @param immuneExposureData immune exposure data.table
#' #' @export
#' #'
#' addImmuneExposureFields <- function(dt, immuneExposure){
#'   pidOrder <-  match(dt$participant_id, expo$participant_id)
#'   dt$exposure_material <- immuneExposure$`Exposure Material Reported`[ pidOrder ]
#'   dt$exposure_process <- immuneExposure$`Exposure Process Preferred`[ pidOrder ]
#'   return(dt)
#' }

#' Update study_time_collected for known issues
#'
#' @param geMetaData meta-data data.table
#' @param studies list of studies to update
#' @param originalTp original timepoint value
#' @param newTp new timepoint value
#' @export
#'
updateStudyTimepoints <- function(geMetaData, studies, originalTp, newTp){
  inTargetStudies <- geMetaData$study_accession %in% studies
  if(originalTp != 24){
    badTimepoint <- geMetaData$time_post_last_vax == originalTp
  }else{
    badTimepoint <- geMetaData$time_post_last_vax == originalTp
  }
  geMetaData$time_post_last_vax[ inTargetStudies & badTimepoint ] <- newTp

  return(geMetaData)
}


#' Remove samples from certain timepoints from desired studies
#'
#' @param esets list of expressionSetObjects
#' @param study study id
#' @param timepoint study_time_collected value
#' @export
#'
removeTimepointFromEset <- function(esets, study, timepoint){
  nms <- names(esets)
  esLoc <- grep(study, names(esets))

  if(timepoint == "Hours"){
    esets[esLoc] <- lapply(esets[esLoc], function(eset){
      eset <- eset[ , eset$study_time_collected_unit != timepoint ]
    })
  }else{
    esets[esLoc] <- lapply(esets[esLoc], function(eset){
      eset <- eset[ , as.numeric(eset$study_time_collected) != timepoint ]
    })
  }

  names(esets) <- nms
  return(esets)
}

#' SDY212 has a biosample that does not have matching meta-data and cannot
#' be explained by the study authors. This function removes that sample.
#'
#' @param esets list of expressionSet objects
#' @export
#'
removeSDY212MissingSample <- function(esets){
  nms <- names(esets)
  esets <- lapply(seq(1:length(esets)), function(index){
    eset <- esets[[index]]
    name <- names(esets)[[index]]
    if(grepl("SDY212", name)){
      eset <- eset[ , eset$biosample_accession != "BS694717"]
    }
    return(eset)
  })
  names(esets) <- nms
  return(esets)
}

#' To accomodate studies with booster shots, add a timepoint post final vaccination
#'
#' @param dt pData from expressionSet object
#' @export
#'
addTimePostLastVax <- function(dt){
  dt$time_post_last_vax <- ifelse(dt$study_accession == "SDY1293",
                                          case_when(
                                            dt$study_time_collected == 60 ~ 0,
                                            dt$study_time_collected == 61 ~ 1,
                                            dt$study_time_collected == 63 ~ 3,
                                            dt$study_time_collected == 74 ~ 14
                                          ),
                                          dt$study_time_collected)
  dt$unit_post_last_vax <- dt$study_time_collected_unit
  return(dt)
}

#' Ensure every subject has baseline gene expression data
#'
#' @param eset expressionSet object
#' @export
#'
removeSubjectsWithoutBaseline <- function(eset){
  pd <- pData(eset)
  allPids <- unique(pd$participant_id)
  pidsWithBaseline <- unique(pd$participant_id[ pd$study_time_collected >= -7 || pd$study_time_collected <= 0 ])
  pidsToRm <- setdiff(allPids, pidsWithBaseline)
  eset <- eset[ , !eset$participant_id %in% pidsToRm ]
}

#' Add matrix, featureset, and cell_type fields to pData objects
#'
#' @param phenoDataSets list of pData objects
#' @param geMatrices ImmuneSpace connection con$cache$GE_matrices object
#' @export
#'
addMatrixRelatedFields <- function(phenoDataSets, geMatrices){
  phenoDataSets <- lapply(seq(1:length(phenoDataSets)), function(index){
    pd <- phenoDataSets[[index]]
    pd$matrix <- geMatrices$name[[index]]
    pd$featureset <- geMatrices$featureset[[index]]
    pd$cell_type <- strsplit(pd$cohort_type, "_")[[1]][[2]]
    return(pd)
  })
}

#' Add name of feature annotation set name used to generate gene expression matrix
#'
#' @param geMetaData single data.table containing all gene expression meta-data
#' @param fasMap feature annotation set map table from ImmuneSpace DB
#' @export
#'
addFeatureAnnotationSetName <- function(geMetaData, featureAnnotationMap){
  geMetaData$featureSetName <- unlist(sapply(geMetaData$featureset, function(x){
    colToUse <- ifelse(x %in% featureAnnotationMap$origid, "origid", "currid")
    return(featureAnnotationMap$name[ featureAnnotationMap[[colToUse]] == x])
  }))
  return(geMetaData)
}

#' Add name of feature annotation set vendor used to generate gene expression matrix
#'
#' @param geMetaData single data.table containing all gene expression meta-data
#' @param fas feature annotation set table from ImmuneSpace DB
#' @export
#'
addFeatureAnnotationSetVendor <- function(geMetaData, featureAnnotation){
  geMetaData$featureSetVendor <- featureAnnotation$vendor[ match(geMetaData$featureSetName, featureAnnotation$name)]
  return(geMetaData)
}

#' Add name of feature annotation set vendor used to generate gene expression matrix
#'
#' @param geMetaData single data.table containing all gene expression meta-data
#' @param gef gene expression files dataset from ImmuneSpace connection object
#' @export
#'
addGSMAccessions <- function(geMetaData, gef){
  gef <- gef[ !is.na(gef$geo_accession), ]
  geMetaData$gsm <- gef$geo_accession[ match(geMetaData$biosample_accession, gef$biosample_accession)]
  return(geMetaData)
}

#' Add demographic variables needed for analysis
#'
#' @param geMetaData single data.table containing all gene expression meta-data
#' @export
#'
addAnalysisVariables <- function(geMetaData){
  geMetaData <- geMetaData[ , c("Hispanic",
                                "study_time_collected",
                                "White",
                                "Asian",
                                "Black",
                                "study_accession2"
                                ) :=
                              list(1 * (ethnicity == 'Hispanic or Latino'),
                                   as.numeric(as.character(study_time_collected)),
                                   1 * (race == 'White'),
                                   1 * (race == 'Asian'),
                                   1 * (race == 'Black or African American'),
                                   ifelse(study_accession == 'SDY1276',
                                          paste(study_accession, gender, sep = '_'),
                                          study_accession))
                            ]
}

#' Subset to just columns needed for cross-study normalization and batch-correction
#'
#' @param geMetaData single data.table containing all gene expression meta-data
#' @export
#'
subsetToOnlyNeededColumns <- function(geMetaData){
  expectedCols <- c("uid",
                    "participant_id",
                    "biosample_accession",
                    "study_time_collected",
                    "study_time_collected_unit",
                    "time_post_last_vax",
                    "unit_post_last_vax",
                    "age_reported",
                    "age_imputed",
                    "gender",
                    "race",
                    "ethnicity",
                    "exposure_material_reported",
                    "exposure_process_preferred",
                    "matrix",
                    "gsm",
                    "study_accession",
                    "study_accession2",
                    "Hispanic",
                    "White",
                    "Asian",
                    "Black",
                    "cell_type",
                    "cohort",
                    "featureSetName",
                    "featureSetVendor"
  )

  geMetaData <- geMetaData[ , ..expectedCols ]
}

#' Summarize probe-level gene expression data by gene symbol selecting the probe
#' for each gene symbol that has maximum average value across all samples within
#' that expression set
#'
#' @param esets list of expressionSet object
#' @export
#'
summarizeByGeneSymbol <- function(esets){
  summarizedEsets <- list()

  for(i in 1:length(esets)){
    # Calc probe average without log transform and DO NOT allow NAs
    # Remove any rows with NAs - SDY1289 montreal cohort 7 probes and SDY212 single probe
    exprs <- data.table(Biobase::exprs(esets[[i]]), keep.rownames = TRUE)
    goodRows <- apply(exprs[,-1], 1, function(row){ all(!is.na(row)) })
    exprs <- exprs[ goodRows ]
    calcAvgWithoutLog <- function(row){ sum(2^row) / length(row) }
    exprs[ , prb_avg := apply(exprs[,-1], 1, calcAvgWithoutLog) ]

    # Latest gene_symbols come from fData based on org.Hs.eg.db - v3.6.0.
    # Ensure that feature data is accurate and ordered correctly before
    # assigning a gene_symbol
    fdat <- Biobase::fData(esets[[i]])
    fdat <- fdat[ fdat$FeatureId %in% exprs$rn, ]
    fdat <- fdat[ order(match(fdat$FeatureId, exprs$rn)), ]
    exprs$gs <- fdat$gene_symbol

    # Check for duplicates at probe level and remove.
    # Gene symbols for these duplicates are noted in IS2_removed_data.
    # This issue is likely caused by gene_symbols being updated in the latest annotation.
    # E.g. probe 1 previously mapped to gene A and gene B, then gene A was updated
    # to be gene B as well.
    exprs <- exprs[ !duplicated(exprs), ]

    # Check for probes mapping to multiple gene symbols and remove
    tmp <- data.table(probe = rownames(exprs),
                      gs = exprs$gs,
                      stringsAsFactors = FALSE)

    probesWithMultipleGS <- tmp[ , cnt := .N, by = "probe"][ cnt > 1 ]

    if(nrow(probesWithMultipleGS) > 0){
      exprs <- exprs[ !(rownames(exprs) %in% probesWithMultipleGS$probe), ]
    }

    # filter to max probe for each gene symbol
    maxPrb <- exprs[, prb_max := max(prb_avg) , by = "gs" ][ prb_avg == prb_max ]

    # Check and remove summary level duplicates
    maxPrb <- maxPrb[ !(duplicated(maxPrb)), ]

    # Check and remove summary level NA values
    maxPrb <- maxPrb[ !is.na(maxPrb$gs) ]

    # handle cases where duplicated gs, but slightly diff vals
    # e.g. SDY1289 - Lausanne Adult Cohort
    if( any(duplicated(maxPrb$gs)) ){
      dupGs <- maxPrb$gs[ duplicated(maxPrb$gs)]
      for(dup in dupGs){
        rows <- grep(dup, maxPrb$gs)
        rmRows <- rows[1:length(rows) - 1]
        maxPrb <- maxPrb[ -rmRows, ]
      }
    }

    maxPrb[ , c("prb_avg", "prb_max", "rn") := NULL ]
    summarizedEsets[[i]] <- maxPrb
  }
  return(summarizedEsets)
}

#' Impute gender based on expression of Y chromosome-related genes
#'
#' @param eset expressionSet
#' @export
#'
imputeGender <- function(eset){

  # Create subset eset with only one sample per participantId.
  # By using !duplicated, we are taking the earliest timepoint
  e0 <- eset[, !duplicated(eset$participant_id) ]
  PD0 <- data.table(pData(e0))
  PD0[, gender_imputed := NA]

  yChromGenes <- intersect(featureNames(e0), yChromGenes)
  yChromGenes <- yChromGenes[ which(rowSums(is.na(exprs(e0)[yChromGenes,])) == 0) ]

  studies <- unique(as.character(e0$study_accession2))
  for(study in studies){
    smpls <- sampleNames(e0)[ which((e0$study_accession2 == study)) ]
    a <- massi_cluster(exprs(e0)[yChromGenes, smpls])
    PD0$gender_imputed[match(smpls, PD0$uid)] <- ifelse(a[[1]]$clustering[smpls] == 1,
                                                        'Female','Male')
  }

  # Check that all PD0 samples have gender_imputed of "Female" or "Male"
  if (!all(PD0$gender_imputed %in% c("Female", "Male"))) {
    stop("gender imputation did not work.")
  }

  # Update eset with imputed gender and age
  PD <- merge(pData(eset),
              subset(PD0, select = c('participant_id','gender_imputed')),
              by = c('participant_id'),
              all.x = TRUE)
  rownames(PD) <- as.character(PD$uid)
  pData(eset) <- PD[sampleNames(eset),]

  return(eset)
}
