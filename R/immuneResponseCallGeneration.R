#' Custom pre-processing for particularly ornery studies
#'
#' @param assay assay name
#' @param df assay data-frame
#' @export
#'
customProcessing <- function(assay, df){

  dt <- data.table(df)
  dt$value_preferred <- as.numeric(dt$value_preferred)

  # Fix SDY1276 log scaling - both assays
  dt <- dt[ dt$study_accession == "SDY1276", value_preferred := 4 ^ dt$value_preferred ]

  if(assay == "neut_ab_titer"){

    # Fix baseline values for SDY1289
    dt <- dt[ dt$study_accession == "SDY1289" &
                dt$value_preferred == 0 &
                as.numeric(dt$study_time_collected) == 0,
              value_preferred := 1]

    # Create baseline data for SDY1264
    sdy1264 <- dt[ dt$study_accession == "SDY1264" ]
    dayZero <- copy(sdy1264)
    dayZero[, dayZero$study_time_collected := '0']
    dayZero[, dayZero$value_preferred := 1 ]
    dayZero[, dayZero$value_reported := '1' ]
    dupes <- which(duplicated(dayZero$participant_id))
    if (length(dupes) > 0) {
      dayZero <- dayZero[ -which ]
    }
    dt <- rbind(dt, dayZero)

  }
  return(dt)
}

#' Pre-processing for 'titer::Calculate' functions
#'
#' @param dt data.table with HAI or NAb data
#' @param postVaxDayRange range of integer values for study_time_collected cutpoints
#' @export
#'
preProcessImmData <- function(dt, postVaxDayRange){

  # remove dupes - coming from ImmuneSpace, not generated
  dt <- dt[ !duplicated(dt) ]

  # Filter data down to samples with selected post-baseline or baseline timepoints
  dt$study_time_collected <- as.numeric(dt$study_time_collected)
  dt <- dt[ study_time_collected %in% seq(postVaxDayRange[[1]], postVaxDayRange[[2]]) |
              study_time_collected <= 0 ]

  # Filter to rows representing max titer value for each timepoint and virus
  dt <- dt[, .SD[value_preferred == max(value_preferred)],
           by = c("participant_id", "virus", "study_time_collected")]

  # define type (post-baseline vs baseline)
  dt$sample_type <- ifelse(dt$study_time_collected <= 0, "pre", "post")

  # Filter baseline samples to Day-0 or closest <0 day
  pre <- dt[ sample_type == "pre",
             .SD[study_time_collected == max(study_time_collected)],
             by = c("participant_id", "virus") ]

  # Filter columns to only those needed and rename as necessary
  pre <- pre[, list(Pre = value_preferred,
                    Study_time_collected_pre = study_time_collected,
                    virus,
                    age_imputed,
                    participant_id,
                    study_accession,
                    irpBatchName)]

  # Filter the post-baseline samples down to those that have the max value_preferred
  # and then if there are ties, use the study_time_collected closest to baseline
  post <- dt[ sample_type == "post",
              .SD[value_preferred == max(value_preferred),
                  .SD[study_time_collected == min(study_time_collected)]
                  ],
              by = c("participant_id", "virus") ]

  # Filter columns to only those needed and rename as necessary
  post <- post[, list(Post = value_preferred,
                      Study_time_collected_post = study_time_collected,
                      virus,
                      age_imputed,
                      participant_id,
                      study_accession,
                      irpBatchName)]

  # Put baseline and post back together - reduces to only those with both pre and post
  full <- merge(pre, post, by = c('participant_id',
                                  'virus',
                                  'study_accession',
                                  'age_imputed',
                                  'irpBatchName'))

  # Update Strain to add study and age_cohort for identifiability
  # full[, study_accession := paste(study_accession,
  #                                 ifelse(full$age_imputed > ageCutoffs[[1]], "old", "young"),
  #                                 sep = "_")]
  full[, virus := paste(irpBatchName, virus, sep = "_")]

  # Ensure that every study * age_cohort * strain has matching numbers of Subjects
  # by removing individuals that do not have values for all strains in their
  # study * age_cohort group
  full <- full[, sumStrains := length(unique(virus)), by = c("irpBatchName")]
  full <- full[, indivStrains := .N, by = c("irpBatchName", "participant_id")]
  full <- full[ sumStrains == indivStrains ]

  # Remove unnecessary columns and update names
  full <- full[, list(SubjectID = participant_id,
                      Strain = virus,
                      study_accession,
                      irpBatchName,
                      Study_time_collected_pre,
                      Pre,
                      Study_time_collected_post,
                      Post)]

  # Split into titer list for FormatTiters based on study
  titer_list_study <- split(full, f = full$irpBatchName)

  # Ensure correct data.frame format for the FormatTiters call
  titer_list_study <- lapply(X = titer_list_study, FUN = as.data.frame)
}



#' Flexible function that wraps around `titer::Calculate` methods while
#' also handling study * age_cohorts that cannot be modelled with an exponential model
#'
#' @param titer_list pre-processed list of dataframes with hai or nab data
#' @param analysis RBA or MFC
#' @param discretize cut points to use for discretizing response call groups
#' @export
#'
performTiterAnalysis <- function(titer_list, analysis, discretize){

  # ----- MAX RESIDUAL BASELINE ADJUSTED -------
  # Discretize using maximum residual after baseline-adjustment for each subject
  # - model each strain using linear model of Pre-vax value ~ Fold Change
  # - get residuals for each subject x strain (create matrix of subs x strain)
  # - get max value of all residuals for each subject
  # - discretize (bin) subjects by the quantile value given (eg. <0.3, 0.3-0.7, >0.7)
  # Calculates the baseline-adjusted fold change for each strain of virus using (unnormalized)
  # fold change and baseline titers. Linear regression or an exponential curve is used to remove
  # the effect of baseline titers on fold changes. The score function (scoreFun) is used to
  # combine the adjusted fold change across multiple strains. Missing (NA) values are handled by
  # being returned as missing in the endpoints in the output
  if( analysis == "RBA" ){

    # Run once one each strain to find strains that do not converge
    # in an exponential model.  It is unclear if there is a programmatic
    # way to determine the conditions that lead to non-convergence prior
    # to running the main function.
    validatedData <- list()
    for(study in names(titer_list)){
      strainsInStudy <- list()
      studyData <- titer_list[[study]]
      for(strain in names(studyData) ){
        ret <- tryCatch(
          titer::Calculate_maxRBA(studyData[strain], discretize = discretize),
          error = function(e) return(e)
        )
        if( "models" %in% names(ret) ){
          validatedData[[study]][[strain]] <- titer_list[[study]][[strain]]
        }else{
          # print(strain)
        }
      }
    }

    # list of non-convergers for reference:

    # ---- HAI ----
    # SDY1119_old_B/Wisconsin/01/2010
    # SDY269_LAIV_young_A/South Dakota/06/2007
    # SDY269_LAIV_young_B/Florida/4/2006
    # SDY400_young_A/California/7/2009

    # ---- NAb ----
    # SDY1264_young_YF17D
    # SDY1289_young_YF17D
    # SDY1294_young_Yellow fever virus 17D
    # SDY180_Pneunomax_young_P. pneumoniae Serotype 12 (12F)
    # SDY180_Pneunomax_young_P. pneumoniae Serotype 4 (4)
    # SDY1325_young_Neisseria meningitidis strain A (F8238)
    # SDY80_old_A/Brisbane/59/2007

    analysisResults <- lapply(validatedData, titer::Calculate_maxRBA, discretize = discretize)

    # ----- MAX FOLD CHANGE -------
  }else{
    # Note that ALL strains, including those that fail RBA are used here
    analysisResults <- lapply(titer_list, titer::Calculate_MFC, discretize = discretize)
  }

  # ------ PREP FOR MERGING INTO RESPONSE DATA SET -----
  # Need to allow for variable discretization values and length(values)
  # Analysis Results come back in form of list of study * age_cohort
  tmp <- lapply(names(analysisResults), function(study){
    studyRes <- analysisResults[[study]]
    analysisTerm <- ifelse(analysis == "RBA",
                           paste0("max", analysis),
                           analysis)
    colsToUse <- sapply(discretize, function(discTerm){
      cname <- paste0(analysisTerm, "_d", gsub("0\\.", "", as.character(discTerm)), "0")
    })

    maxStrainTitle <- paste0("maxStrain_", analysis)
    baselineTimepoint <- paste0("ImmResp_baseline_timepoint_", analysis)
    baselineValue <- paste0("ImmResp_baseline_value_", analysis)
    postVaxTimepoint <- paste0("ImmResp_postVax_timepoint_", analysis)
    postVaxValue <- paste0("ImmResp_postVax_value_", analysis)

    resultsList <- list()

    resultsList$SubjectID <- names(studyRes[[analysisTerm]])
    resultsList[[analysisTerm]] <- studyRes[[analysisTerm]]

    newCols <- c(baselineTimepoint, baselineValue, postVaxTimepoint, postVaxValue)
    for(cname in newCols){
      resultsList[[ cname ]] <- rep(NA, length(resultsList$SubjectID))
    }

    for(colNm in colsToUse){
      resultsList[[gsub("d", "p", colNm)]] <- studyRes[[colNm]]
    }

    if( analysis == "RBA"){
      if(length(colnames(studyRes$residualMatrix)) == 1){
        resultsList[[maxStrainTitle]] <- rep(colnames(studyRes$residualMatrix)[[1]],
                                             length(resultsList$SubjectID))
      }else{
        resultsList[[maxStrainTitle]] <- unlist(sapply(seq(1, length(studyRes$maxRBA)), function(i){
          maxValue <- studyRes$maxRBA[[i]]
          strain <- names(which(studyRes$residualMatrix[i, ] == maxValue))
        }))
      }
    }else{
      strainFCs <- data.frame(sapply(titer_list[[study]], '[', 'FC'))
      colnames(strainFCs) <- names(titer_list[[study]])
      resultsList[[maxStrainTitle]] <- apply(strainFCs, 1, function(i){
        strain <- names(which.max(i))
      })
    }

    for(i in seq(1, length(resultsList$SubjectID))){
      strainInfoDF <- titer_list[[study]][[ resultsList[[maxStrainTitle]][[i]] ]]
      rowId <- which(strainInfoDF$SubjectID == resultsList$SubjectID[[i]])
      resultsList[[baselineTimepoint]][[i]] <- strainInfoDF$Study_time_collected_pre[ rowId ]
      resultsList[[postVaxTimepoint]][[i]] <- strainInfoDF$Study_time_collected_post[ rowId ]
      resultsList[[baselineValue]][[i]] <- strainInfoDF$Pre[ rowId ]
      resultsList[[postVaxValue]][[i]] <- strainInfoDF$Post[ rowId ]
    }

    resultsList[[maxStrainTitle]] <- gsub("SDY\\d{2,4}.*_", "", resultsList[[maxStrainTitle]])

    ret <- data.frame(resultsList, stringsAsFactors = FALSE)
  })

  tmp <- rbindlist(tmp)
}

#' run both maxRBA and MFC calculations using the 'titer' package by
#' Stefan Avey of Yale University
#'
#' @param titer_list pre-processed list of dataframes with hai or nab data
#' @param df immdata df
#' @param discretizationValues cut points to use for discretizing response call groups
#' @export
#'
runAllAnalyses <- function(titer_list, df, discretizationValues){

  # titer:: Calculate_maxRBA
  tmp_rba <- performTiterAnalysis(titer_list = titer_list,
                                  analysis = "RBA",
                                  discretize = discretizationValues[["RBA"]])
  rba_df <- merge(df, tmp_rba,
              by.x = "participant_id", by.y = 'SubjectID')

  # titer::Calculate_MFC (multiple fold change)
  tmp_mfc <- performTiterAnalysis(titer_list = titer_list,
                                  analysis = "MFC",
                                  discretize = discretizationValues[["MFC"]])
  mfc_df <- merge(df, tmp_mfc,
              by.x = "participant_id", by.y = 'SubjectID')

  sharedCols <- intersect(colnames(rba_df), colnames(mfc_df))
  all <- merge(mfc_df, rba_df, by = sharedCols, all = TRUE) # ensure those with MFC, not RBA stick around
  all <- all[ !duplicated(all) ]
}

#' Generate immune response calls for HAI or NAb using pipeline originally
#' developed by Daniel Chawla at Yale University
#'
#' @param assay assay name
#' @param df immdata df
#' @param postVaxDayRange Allowable timepoints for post-vaccine values
#' @param discretizationValues cut points to use for discretizing response call groups
#' @export
#'
generateNAbHAIresponse <- function(assay, df, postVaxDayRange, discretizationValues){
  dt <- customProcessing(assay = assay,
                         df = df)

  titer_list_study <- preProcessImmData(dt = dt,
                                        postVaxDayRange = postVaxDayRange)

  titer_list <- suppressMessages(lapply(X = titer_list_study,
                                        FUN = titer::FormatTiters,
                                        log2Transform = TRUE,
                                        fcMinZero = FALSE))

  df <- runAllAnalyses(titer_list = titer_list,
                       df = df,
                       discretizationValues = discretizationValues)
}

#' Generate immune response call for age specific ELISA cohort
#'
#' @param dt immdata df
#' @param discretizationValues points to use for cutting of low, mid, high responders
#' @param postVaxDayRange range of allowable study_time_collected
#' @export
#'
generateELISAResponse <- function(dt, discretizationValues, postVaxDayRange){

  dt <- dt[ grep("IgG|^Hepatitis", dt$analyte) ]

  # Filter data down to samples with selected post-baseline or baseline timepoints
  dt$study_time_collected <- as.numeric(dt$study_time_collected)

  # Correct SDY1328 dates per comments of S. Fourati - post-vax in month 2
  datesToChange <- dt$study_accession == "SDY1328" & dt$study_time_collected == 7
  dt$study_time_collected[ datesToChange] <- 30

  # Subset
  postVaxTp <- seq(postVaxDayRange[[1]], postVaxDayRange[[2]])
  dt <- dt[ dt$study_time_collected %in% c(0, postVaxTp) ]
  dt$value_preferred <- as.numeric(dt$value_preferred)

  # SDY1260 Corrections
  dt$analyte[ grep("IgG(\\d|)_serotype", dt$analyte)] <- "IgG"
  dt$value_preferred[dt$study_accession == "SDY1260"] <- `^`(2, dt$value_preferred[dt$study_accession == "SDY1260"])

  # SDY1328 - ensure day 0 are considered "naive"
  # Standardize the FC for naive subjects that show no change to be 0
  samplesToUpdate <- dt$study_accession == "SDY1328" &
    ( dt$value_preferred == 2.5 | dt$study_time_collected == 0)
  dt$value_preferred[ samplesToUpdate ] <- 1

  # Only applies to SDY1260 - sum Serotype A and Serotype C
  dt <- dt[, dt$value_preferred := sum(value_preferred),
           by = c("participant_id", "study_time_collected", "vaccine", "vaccine_type", "pathogen")]

  colsCreatingDupes <- c("expsample_accession", "value_reported", "unit_reported")
  dt <- dt[, c(colsCreatingDupes) := NULL ]

  # De-dupe for SDY1260
  dt <- dt[!duplicated(dt)]

  # define type (post-baseline vs baseline)
  dt$sample_type <- ifelse(dt$study_time_collected == 0, "pre", "post")

  # Filter baseline samples to Day-0 or closest <0 day
  pre <- dt[ dt$sample_type == "pre"]

  # Filter columns to only those needed and rename as necessary
  pre <- pre[, c("ImmResp_baseline_value_MFC",
                 "ImmResp_baseline_timepoint_MFC",
                 "maxStrain_MFC")
             :=
               list(value_preferred,
                    study_time_collected,
                    analyte)
             ]

  # Filter the post-baseline samples down to those that have the max value_preferred
  post <- dt[ sample_type == "post"]

  # Filter columns to only those needed and rename as necessary
  post <- post[, c("ImmResp_postVax_value_MFC",
                   "ImmResp_postVax_timepoint_MFC",
                   "maxStrain_MFC")
               :=
                 list(value_preferred,
                      study_time_collected,
                      analyte)
               ]

  colsToRm <- c("biosample_accession",
                "study_time_collected",
                "value_preferred",
                "uid",
                "sample_type",
                "comments",
                "description",
                "phenotype",
                "unit_preferred")

  pre <- pre[, c(colsToRm) := NULL ]
  post <- post[, c(colsToRm) := NULL ]

  # Put baseline and post back together - reduces to only those with both pre and post
  sharedCols <- intersect(colnames(post), colnames(pre))
  full <- merge(pre, post, by = sharedCols, all = TRUE)

  # only those with both pre and post
  full <- full[ !is.na(full$ImmResp_postVax_value_MFC) & !is.na(full$ImmResp_baseline_value_MFC) ]

  # add MFC - Actually using baseline as background and subtracting instead of doing proper MFC
  full$MFC <- as.numeric(full$ImmResp_postVax_value_MFC) - as.numeric(full$ImmResp_baseline_value_MFC)

  # discretize within study
  discretize <- function(values, cutPoint){
    x <- stats::quantile(values, c(cutPoint, 1 - cutPoint))
    res <- sapply(values, function(y){
      if(y <= x[[1]]){
        return("lowResponder")
      }else if(y >= x[[2]]){
        return("highResponder")
      }else{
        return("moderateResponder")
      }
    })
  }

  for(point in discretizationValues){
    colName <- paste0("MFC_p", gsub("0\\.", "", point), "0")
    full[ , c(colName) := discretize(MFC, point), by = "study_accession"]
  }

  full[, analyte := NULL ]

  return(full)
}

#' Generate immune response calls for HAI, NAb or ELISA
#'
#' @param assay assay name
#' @param data immdata df
#' @param postVaxDayRange Allowable timepoints for post-vaccine values
#' @param discretizationValues cut points to use for discretizing response call groups
#' @export
#'
generateResponseCall <- function(assay, data, postVaxDayRange, discretizationValues){
  if(assay %in% c("hai","neut_ab_titer")){
    res <- generateNAbHAIresponse(
      assay = assay,
      df = data,
      postVaxDayRange = postVaxDayRange,
      discretizationValues = discretizationValues
    )
  }else if(assay == "elisa"){
    res <- generateELISAResponse(
      dt = data,
      discretizationValues = discretizationValues$MFC,
      postVaxDayRange = postVaxDayRange
    )
  }
  return(res)
}
