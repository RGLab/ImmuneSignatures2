#' Add fields from immune exposure data to meta-data
#'
#' @param dt meta-data data.table
#' @param immuneExposureData immune exposure data.table
#' @export
#'
addImmuneExposureFields <- function(dt, immuneExposure){
  pidOrder <-  match(dt$participant_id, expo$participant_id)
  dt$exposure_material <- immuneExposure$`Exposure Material Reported`[ pidOrder ]
  dt$exposure_process <- immuneExposure$`Exposure Process Preferred`[ pidOrder ]
  return(dt)
}

#' Update study_time_collected for known issues
#'
#' @param geMetaData meta-data data.table
#' @param studies list of studies to update
#' @param originalTp original timepoint value
#' @param newTp new timepoint value
#' @export
#'
updateStudyTimepoints <- function(geMetaData, studies, originalTp, newTp){
  if(originalTp != 24){
    geMetaData$study_time_collected[ geMetaData$study %in% studies,
                                     geMetaData$study_time_collected == originalTp] <- newTp
  }else{
    geMetaData$study_time_collected[ geMetaData$study %in% studies,
                                     geMetaData$study_time_collected >= originalTp] <- newTp
  }

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

  return(esets)
}

summarizeByGeneSymbol <- function(esets){
  summarizedEsets <- list()

  for(i in 1:length(esets)){
    # Calc probe average without log transform and DO NOT allow NAs
    # Remove any rows with NAs - SDY1289 montreal cohort 7 probes and SDY212 single probe
    exprs <- data.frame(Biobase::exprs(esets[[i]]))
    goodRows <- apply(exprs, 1, function(row){ all(!is.na(row)) })
    exprs <- exprs[ goodRows, ]
    exprs$prb_avg <- apply(exprs, 1 , function(row){ sum(2^row) / length(row) })

    # Latest gene_symbols come from fData based on org.Hs.eg.db - v3.6.0.
    # Ensure that feature data is accurate and ordered correctly before
    # assigning a gene_symbol
    fdat <- Biobase::fData(esets[[i]])
    fdat <- fdat[ fdat$FeatureId %in% rownames(exprs), ]
    fdat <- fdat[ order(match(fdat$FeatureId, rownames(exprs))), ]
    exprs$gs <- fdat$gene_symbol

    # Check for duplicates at probe level and remove.
    # Gene symbols for these duplicates are noted in IS2_removed_data.
    # This issue is likely caused by gene_symbols being updated in the latest annotation.
    # E.g. probe 1 previously mapped to gene A and gene B, then gene A was updated
    # to be gene B as well.
    exprs <- exprs[ !duplicated(exprs), ]

    # Check for probes mapping to multiple gene symbols and remove
    tmp <- data.frame(probe = rownames(exprs),
                      gs = exprs$gs,
                      stringsAsFactors = FALSE)

    probesWithMultipleGS <- tmp %>%
      group_by(probe) %>%
      summarize(cnt = n()) %>%
      filter(cnt > 1)

    exprs <- exprs[ !(rownames(exprs) %in% probesWithMultipleGS$probe), ]

    # filter to max probe for each gene symbol
    maxPrb <- exprs %>%
      group_by(gs) %>%
      filter(prb_avg == max(prb_avg)) %>%
      ungroup()

    # Check and remove summary level duplicates
    maxPrb <- maxPrb[ !(duplicated(maxPrb)), ]

    # Check and remove summary level NA values
    maxPrb <- data.frame(maxPrb[ !is.na(maxPrb$gs), ])

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

    summarizedEsets[[i]] <- maxPrb[, colnames(maxPrb) != "prb_avg" ]
  }
  return(summarizedEsets)
}
