#' Generate sample MDS plots for QC
#'
#' @param eset expressionSet
#' @export
#'
qualityControl.sampleMDSPlot <- function(eset, numberOfSamples = 10, colorCol = "study_accession"){
  # Based on code by Aris - 03/10/2019
  day0 <- eset[ , eset$time_post_last_vax == 0 ]
  studies <- unique(day0$study_accession)

  ind <- unlist(lapply(studies, function(study){
    smpls <- sampleNames(day0)[ day0$study_accession == study ]
    if(length(smpls) > numberOfSamples){
      smpls <- sample(smpls, numberOfSamples)
    }
    return(smpls)
  }))

  subsetDay0 <- day0[ , sampleNames(day0) %in% ind ]
  colors <- RColorBrewer::brewer.pal(10, "Spectral")

  colorVec <- unique(subsetDay0[[colorCol]])
  colors <- colorRampPalette(colors)(length(colorVec))
  tmp <- plotMDS(subsetDay0, col = colors[factor(colorVec)], labels = colorVec)
}

#' Generate sample MDS plots for QC
#'
#' @param eset expressionSet
#' @param returnObject options are allMatricesPlot (default), probSamplesPlot, probSamplesDT
#' @import ggbeeswarm tidyverse Biobase data.table
#' @export
#'
qualityControl.genderByMatrix <- function(eset, returnObject = "allMatricesPlot"){
  yChromGenes <- yChromGenes[ yChromGenes %in% featureNames(eset) ]

  plotDF <- colMeans(exprs(eset)[yChromGenes, ], na.rm = TRUE) %>%
    data.frame(chry = .) %>%
    rownames_to_column() %>%
    merge(y    = pData(eset),
          by.x = "rowname",
          by.y = "uid")

  if(returnObject == "probSubjects"){
    colors <- RColorBrewer::brewer.pal(12, "Paired")
    set.seed(10)
    colors <- sample(colors, length(colors))
    colors <- colorRampPalette(colors)(length(unique(plotDF$participant_id)))

    fullPlot <- ggplot(data = plotDF,
                       mapping = aes(x = gender_imputed, y = chry)) +
      geom_boxplot(outlier.color = "transparent", fill = "grey") +
      geom_jitter(height = 0, width = 0.25, mapping = aes(color = participant_id)) +
      labs(y = "Average probe intensities (chrY)") +
      scale_color_manual(values = colors) +
      facet_wrap(facets = ~study_accession+matrix, scale = "free") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.pos = "none",
            strip.text = element_text(size = 5))
    return(fullPlot)
  }

  # find outlier (1.5 x IGR from Q1 and Q3)
  outlierDF <- plotDF %>%
    group_by(study_accession, matrix, gender_imputed) %>%
    mutate(up = chry >= quantile(chry, probs = 0.75) + 1.5 * IQR(chry),
           dn = chry <= quantile(chry, probs = 0.25) - 1.5 * IQR(chry))

  quartileDF <- plotDF %>%
    group_by(study_accession, matrix, gender_imputed) %>%
    summarize(q1 = quantile(chry, probs = 0.25) - 1.5 * IQR(chry),
              q3 = quantile(chry, probs = 0.75) + 1.5 * IQR(chry)) %>%
    mutate(gender_imputed = c("Female" = "Male", "Male" = "Female")[gender_imputed])

  # flag outlier samples (possible swap)
  flagDF <- merge(x  = outlierDF,
                  y  = quartileDF,
                  by = c("study_accession", "matrix", "gender_imputed")) %>%
    mutate(flag = FALSE,
           flag = ifelse(test = gender_imputed %in% "Female" & up & chry >= q1,
                         yes  = TRUE,
                         no   = flag),
           flag = ifelse(test = gender_imputed %in% "Male" & dn & chry <= q3,
                         yes  = TRUE,
                         no   = flag))

  plotDF <- merge(x     = plotDF,
                  y     = select(flagDF, rowname, flag),
                  by    = "rowname",
                  all.x = TRUE) %>%
    # if gender_imputed not specified, set flag as false
    mutate(flag = ifelse(test = is.na(flag),
                         yes  = FALSE,
                         no   = flag))



  if(returnObject == "allMatricesPlot"){
    fullPlot <- ggplot(data = plotDF,
                       mapping = aes(x = gender_imputed, y = chry)) +
      geom_boxplot(outlier.color = "transparent", fill = "grey") +
      geom_jitter(height = 0, width = 0.25, mapping = aes(color = flag)) +
      labs(y = "Average probe intensities (chrY)") +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      facet_wrap(facets = ~study_accession+matrix, scale = "free") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.pos = "none",
            strip.text = element_text(size = 5))
    return(fullPlot)

  }else if(returnObject == "probSamplesDT"){
    return(plotDF[ plotDF$flag, ])

  }else if(returnObject == "probSamplesPlot"){
    problemMatrix <- unique(plotDF$matrix[ plotDF$flag])
    plotTemp <- filter(plotDF, matrix %in% problemMatrix)

    prbSmplsPlot <- ggplot(data = plotTemp,
                           mapping = aes(x = gender_imputed, y = chry)) +
      geom_boxplot(outlier.color = "transparent", fill = "grey") +
      geom_beeswarm(mapping = aes(color = flag), cex = 2.5, size = 0.8) +
      labs(y = "Average probe intensities (chrY)", x = "Sex") +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      facet_wrap(facets = ~study_accession+matrix, scale = "free") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.pos = "none",
            strip.text = element_text(size = 6))
    return(prbSmplsPlot)
  }else{
    stop("returnObject parameter value not recognized")

  }
}

#' Generate table of studies with count of participants from vector of pids
#'
#' @param participantIds
#' @export
#'
qualityControl.createSubjectsByStudyTable <- function(participantIds){
  uniquePids <- unique(participantIds)
  studies <- gsub("SUB\\d{6}\\.", "SDY", uniquePids)
  studiesTbl <- table(studies)
}
