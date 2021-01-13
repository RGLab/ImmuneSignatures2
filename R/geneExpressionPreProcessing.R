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
  pidsWithBaseline <- unique(pd$participant_id[ pd$time_post_last_vax >= -7 & pd$time_post_last_vax <= 0 ])
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
  phenoDataSets <- lapply(names(phenoDataSets), function(matrixName){
    pd <- phenoDataSets[[matrixName]]
    pd$matrix <- matrixName
    pd$featureset <- geMatrices$featureset[[match(matrixName, geMatrices$name)]]
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

# TODO: should this change to complete.cases?
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
    exprs <- data.table(Biobase::exprs(esets[[i]]), keep.rownames = TRUE)
    goodRows <- apply(exprs[,-1], 1, function(row){ all(!is.na(row)) })
    exprs <- exprs[ goodRows ]
    calcAvgWithoutLog <- function(row){ sum(2^row) / length(row) }
    exprs[ , prb_avg := apply(exprs[,-1], 1, calcAvgWithoutLog) ]

    # Ensure that feature data is accurate and ordered correctly before
    # assigning a gene_symbol
    fdat <- Biobase::fData(esets[[i]])
    fdat <- fdat[ fdat$FeatureId %in% exprs$rn, ]
    fdat <- fdat[ order(match(fdat$FeatureId, exprs$rn)), ]
    exprs$gs <- fdat$gene_symbol

    # Check for duplicates at probe level and remove.
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

    # Check and remove summary level duplicates w
    # here multiple probes have exact same average value AND same sample values
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

#' Impute presence of y chromosome based on expression of Y chromosome-related genes
#' and flag for QC if values do not match
#'
#' @param eset expressionSet
#' @export
#'
imputeYchrom.useAllTimepoints <- function(eset){
  PD <- data.table(pData(eset))
  PD[, y_chrom_present_timepoint := NA]

  yGenes <- intersect(featureNames(eset), yChromGenes)
  yGenes <- yGenes[ which(rowSums(is.na(exprs(eset)[yGenes,])) == 0) ]

  matrices <- unique(as.character(eset$matrix))
  for(matrix in matrices){
    print(matrix)

    # SDY1119 TD2 OLD - only one pid, SDY1276 cohorts are separated by gender
    if(!grepl("SDY1276", matrix)){
      smpls <- which(eset$matrix == matrix)
      if(length(smpls) > 2){
        yChromEset <- eset[yGenes, smpls]

        # Get clusters from reduced dimensional projection to eliminate
        # issues
        mds <- limma::plotMDS(yChromEset)
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
          y_chrom_present <- ifelse(mdsCall == 1, TRUE, FALSE)
        }else{
          y_chrom_present <- ifelse(mdsCall == 1, FALSE, TRUE)
        }

        PD$y_chrom_present_timepoint[ match(names(y_chrom_present), PD$uid) ] <- y_chrom_present
      }
    }
  }

  # If ychrom is not imputed based on expression level, impute from reported gender
  ychromImputationNotDone <- which(is.na(PD$y_chrom_present_timepoint))
  PD$y_chrom_present_timepoint[ychromImputationNotDone] <-
    ifelse(PD$gender[ychromImputationNotDone] == "Male", TRUE, FALSE)

  # if y_chrom_present in agreement for all timepoints, use that
  # If y_chrom_present not in agreement for all timepoints:
  # 1. Has gender_reported and all gender_reported in agreement: derive from gender_reported
  # 2. No gender_reported, use majority y_chrom_present

  summarizeYchrom <- function(gender_reported, y_chrom_present_timepoint){
    reported_gender <- unique(gender_reported)
    y_chrom_present <- unique(y_chrom_present_timepoint)

    if (length(y_chrom_present) == 1) {
      # y_chrom_present in agreement for all timepoints
      return(y_chrom_present)
    } else {
      # y_chrom_present not in agreement for all timepoints
      if (length(reported_gender) == 1 & reported_gender %in% c("Male", "Female")) {
        # Has gender_reported an all in agreement
        y_chrom_present <- ifelse(reported_gender == "Male", TRUE, FALSE)
        return(y_chrom_present)
      } else {
        # No gender_reported or not in agreement
        res <- table(y_chrom_present_timepoint)
        majority <- names(res)[ res == max(res) ]
        if(length(majority) == 1){
          # Clear majority
          return(majority)
        }else{
          # No clear majority
          return(y_chrom_present_timepoint)
        }
      }
    }
  }

  # flag problem samples with following logic:
  # 1. All values for y_chrom_present are in agreement: all pass
  # 2. y_chrom_present varies by timepoint:
  #   2a. has gender_reported: derive ychrom from gender_reported and pass if in agreement
  #   2b. no gender_reported and clear majority for ychrom: pass if in agreement with majority
  #   2c. no gender_reported and no clear majority for ychrom: all fail.
  flagProblemTimepoints <- function(gender_reported, y_chrom_present_timepoint){
    if (length(unique(y_chrom_present_timepoint)) == 1) {
      return(FALSE)
    } else if (unique(gender_reported) %in% c("Female", "Male")) {
      # Imputed ychrom values differ by timepoint: Flag when ychrom does not align with
      # expected ychrom value based on gender
      expected_y_chrom <- ifelse(gender_reported == "Male", TRUE, FALSE)
      return(y_chrom_present_timepoint != expected_y_chrom)
    } else {
      # ychrom varies by timepoint and gender_reported not reported:
      # Check if clear majority of ychrom values
      res <- table(y_chrom_present_timepoint)
      majority <- names(res)[ res == max(res) ]
      if (length(majority) == 1) {
        # Clear majority: fail QC if disagrees with majority
        return(y_chrom_present_timepoint != majority)
      } else {
        # No clear majority: all fail.
        return(TRUE)
      }
    }
  }



  PD[ , y_chrom_present := summarizeYchrom(gender, y_chrom_present_timepoint), by = "participant_id"]

  PD[ , failedYchromQC := flagProblemTimepoints(gender, y_chrom_present_timepoint), by = "participant_id"]

  pData(eset) <- PD

  return(eset)
}

