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
    badTimepoints <- geMetaData$time_post_last_vax == originalTp
  }else{
    badTimepoints <- geMetaData$time_post_last_vax >= originalTp
  }
  geMetaData$time_post_last_vax[ inTargetStudies & badTimepoints ] <- newTp

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
  pidsWithBaseline <- unique(pd$participant_id[ pd$study_time_collected >= -7 & pd$study_time_collected <= 0 ])
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
    pd$cell_type <- regmatches(pd$cohort_type,
                               regexpr("Whole blood|PBMC", pd$cohort_type, ignore.case = TRUE))
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
  geMetaData$featureSetVendor <- featureAnnotation$vendor[ match(geMetaData$featureSetName,
                                                                 featureAnnotation$name)]
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
  geMetaData$gsm <- gef$geo_accession[ match(geMetaData$biosample_accession,
                                             gef$biosample_accession)]
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
                                "Black"
                                ) :=
                              list(1 * (ethnicity == 'Hispanic or Latino'),
                                   as.numeric(as.character(study_time_collected)),
                                   1 * (race == 'White'),
                                   1 * (race == 'Asian'),
                                   1 * (race == 'Black or African American'))
                            ]
}

#' Create new featureSetName field that coalesces similar platforms for analysis
#'
#' @param dt data table
#' @export
#'
addCoalescedFeatureSetName <- function(dt){
  dt$featureSetName2 <- dt$featureSetName
  dt$featureSetName2[ grepl('ustom', dt$featureSetName) ] <- 'RNA-seq'
  dt$featureSetName2[ grepl('HT-12', dt$featureSetName) ] <- 'HumanHT-12_2018'
  return(dt)
}

#' Remove incomplete rows
#'
#' @param eset expressionSet
#' @export
#'
removeAllNArows <- function(eset){
  em <- Biobase::exprs(eset)
  allNARows <- apply(em, 1, function(x){ all(is.na(x)) })
  eset <- eset[ !allNARows, ]
}

#' Subset to just columns needed for cross-study normalization and batch-correction
#'
#' @param geMetaData single data.table containing all gene expression meta-data
#' @export
#'
subsetToOnlyNeededColumns <- function(geMetaData){
  geMetaData <- geMetaData[ , ..expectedGeMetaDataColumns ]
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

    # handle cases where duplicated gs due to same average value,
    # but slightly different sample values
    # e.g. SDY1289 - Lausanne Adult Cohort
    if( any(duplicated(maxPrb$gs)) ){
      dupGs <- maxPrb$gs[ duplicated(maxPrb$gs) ]
      for(dup in unique(dupGs)){
        rows <- grep(dup, maxPrb$gs)
        rmRows <- rows[ 1:length(rows) - 1 ]
        maxPrb <- maxPrb[ -rmRows ]
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
imputeGender.useBaseline <- function(eset){

  # Create subset eset with only one sample per participantId, using d0 timepoint preferably
  # Removing technical replicates (only one)
  e0 <- eset[, eset$time_post_last_vax <= 0 ]
  pd <- data.table(pData(e0))
  pdBaseline <- pd[ , maxBaseline := max(time_post_last_vax) , by = "participant_id" ]
  pdBaseline <- pdBaseline[ time_post_last_vax == maxBaseline ]
  pdBaseline <- pdBaseline[ !duplicated(participant_id) ]
  e0 <- e0[, e0$uid %in% pdBaseline$uid ]

  PD0 <- data.table(pData(e0))
  PD0[, gender_imputed := NA]

  yChromGenes <- intersect(featureNames(e0), yChromGenes)
  yChromGenes <- yChromGenes[ which(rowSums(is.na(exprs(e0)[yChromGenes,])) == 0) ]

  # SDY1119 TD2 OLD - only one pid, SDY1276 cohorts are separated by gender
  matrices <- unique(as.character(e0$matrix))
  for(matrix in matrices){
    if(!grepl("SDY1276", matrix)){
      smpls <- sampleNames(e0)[ which((e0$matrix == matrix)) ]
      if(length(smpls) >= 2){
        yChromEset <- e0[yChromGenes, smpls]

        # Get clusters from reduced dimensional projection to eliminate
        # issues
        mds <- plotMDS(yChromEset)
        mdsClust <- kmeans(mds$x, 2)
        mdsCall <- mdsClust$cluster

        # Figure out which of original clusters has higher overall values
        # for expression of yChromGenes
        clustANames <- names(mdsCall)[ mdsCall == 1 ]

        clustAexprs <- yChromEset[ , yChromEset$uid %in% clustANames]
        clustAMeans <- mean(rowMeans(exprs(clustAexprs)))

        clustBexprs <- yChromEset[ , !yChromEset$uid %in% clustANames]
        clustBMeans <- mean(rowMeans(exprs(clustBexprs)))

        if(clustAMeans > clustBMeans){
          mdsFinal <- ifelse(mdsCall == 1, "Male", "Female")
        }else{
          mdsFinal <- ifelse(mdsCall == 1, "Female", "Male")
        }

        PD0$gender_imputed[match(names(mdsFinal), PD0$uid)] <- mdsFinal
      }
    }
  }

  genderImputationFailed <- which(is.na(PD0$gender_imputed))
  PD0$gender_imputed[genderImputationFailed] <- PD0$gender[genderImputationFailed]

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
  orderMatch <- order( match(rownames(PD), colnames(exprs(eset))) )
  pData(eset) <- PD[ orderMatch, ]

  return(eset)
}

#' Impute gender based on expression of Y chromosome-related genes
#'
#' @param eset expressionSet
#' @export
#'
imputeGender.useAllTimepoints <- function(eset){
  PD <- data.table(pData(eset))
  PD[, gender_imputed_timepoint := NA]

  yGenes <- intersect(featureNames(eset), yChromGenes)
  yGenes <- yGenes[ which(rowSums(is.na(exprs(eset)[yGenes,])) == 0) ]

  # SDY1119 TD2 OLD - only one pid, SDY1276 cohorts are separated by gender
  matrices <- unique(as.character(eset$matrix))
  for(matrix in matrices){
    print(matrix)
    if(!grepl("SDY1276", matrix)){
      smpls <- sampleNames(eset)[ which((eset$matrix == matrix)) ]
      if(length(smpls) > 2){
        yChromEset <- eset[yGenes, smpls]

        # Get clusters from reduced dimensional projection to eliminate
        # issues
        mds <- plotMDS(yChromEset)
        mdsClust <- kmeans(mds$x, 2)
        mdsCall <- mdsClust$cluster

        # Figure out which of original clusters has higher overall values
        # for expression of yChromGenes
        clustANames <- names(mdsCall)[ mdsCall == 1 ]

        clustAexprs <- yChromEset[ , yChromEset$uid %in% clustANames]
        clustAMeans <- mean(colMeans(exprs(clustAexprs)))

        clustBexprs <- yChromEset[ , !yChromEset$uid %in% clustANames]
        clustBMeans <- mean(colMeans(exprs(clustBexprs)))

        if(clustAMeans > clustBMeans){
          mdsFinal <- ifelse(mdsCall == 1, "Male", "Female")
        }else{
          mdsFinal <- ifelse(mdsCall == 1, "Female", "Male")
        }

        PD$gender_imputed_timepoint[ match(names(mdsFinal), PD$uid) ] <- mdsFinal
      }
    }
  }


  genderImputationNotDone <- which(is.na(PD$gender_imputed_timepoint))
  PD$gender_imputed_timepoint[genderImputationNotDone] <- PD$gender[genderImputationNotDone]

  # Only reassign if all timepoints are opposite sex
  summarizeGender <- function(genderReported, genderGuesses){
    originalGender <- unique(genderReported)
    assignedGender <- unique(genderGuesses)
    if(length(assignedGender) == 1){
      return(assignedGender)
    }else{
      return(originalGender)
    }
  }

  # flag problem samples if there is disagreement within subject
  flagProblemTimepoints <- function(gender_imputed_timepoint){
    status <- length(unique(gender_imputed_timepoint)) > 1
  }

  PD[ , gender_imputed := summarizeGender(gender, gender_imputed_timepoint), by = "participant_id"]

  PD[ , failedGenderQC := flagProblemTimepoints(gender_imputed_timepoint), by = "participant_id"]

  pData(eset) <- PD

  return(eset)
}

#' Remove selected participants from expression set
#'
#' @param eset expressionSet
#' @param participantIdsToRm vector of participant IDs to remove
#' @export
#'
adjustProblemStudies <- function(eset, problemStudies){
  eset$failedGenderQC[ eset$study_accession %in% problemStudies] <- TRUE
  return(eset)
}

#' Remove selected participants from expression set
#'
#' @param eset expressionSet
#' @param participantIdsToRm vector of participant IDs to remove
#' @export
#'
removeSubjectsWithoutGenderConsensus <- function(eset, participantIdsToRm){
  eset <- eset[ , !eset$participant_id %in% participantIdsToRm]
}
