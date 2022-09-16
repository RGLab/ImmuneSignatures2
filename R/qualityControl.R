#' Generate sample MDS plots for QC
#'
#' @param eset expressionSet
#' @param method MDS or PCA
#' @param numberOfSamples number of samples to use per study
#' @param colorCol field name to use for labeling samples in plot
#' @export
#'
qualityControl.samplePlot <- function(eset,
                                      method = "MDS",
                                      numberOfSamples = 10,
                                      colorCol = "study_accession"){
  # Based on code by Aris - 03/10/2019
  day0 <- eset[ , eset$time_post_last_vax <= 0 ]
  studies <- unique(day0$study_accession)

  ind <- unlist(lapply(studies, function(study){
    smpls <- sampleNames(day0)[ day0$study_accession == study ]
    if(length(smpls) > numberOfSamples){
      smpls <- sample(smpls, numberOfSamples)
    }
    return(smpls)
  }))

  subsetDay0 <- day0[ , sampleNames(day0) %in% ind ]

  if(method == "MDS"){
    colors <- RColorBrewer::brewer.pal(10, "Spectral")
    colorVec <- unique(subsetDay0[[colorCol]])
    colors <- sample(colors, length(colorVec), replace = TRUE)
    colors <- colors[ match(subsetDay0[[colorCol]], colorVec)]
    tmp <- limma::plotMDS(subsetDay0, col = colors, labels = subsetDay0[[colorCol]])
    return(tmp)
  }else if(method == "PCA"){
    em <- as.matrix(exprs(subsetDay0))
    em <- em[ stats::complete.cases(em), ]
    tem <- data.frame(t(em), stringsAsFactors = FALSE)
    res <- stats::prcomp(tem)
    pd <- pData(subsetDay0)
    ggplot2::autoplot(res, data = pd, colour = colorCol)
  }
}

#' Generate box plots showing ychrom imputed after QC
#'
#' @param eset expressionSet
#' @param returnObject options are allMatricesPlot (default), probSamplesPlot, probSamplesDT
#' @export
#'
qualityControl.yChromPresentByMatrix <- function(eset, returnObject = "allMatricesPlot"){
  yChromGenes <- yChromGenes[ yChromGenes %in% featureNames(eset) ]

  plotDF <- colMeans(exprs(eset)[yChromGenes, ], na.rm = TRUE) %>%
    data.frame(chry = .) %>%
    tibble::rownames_to_column() %>%
    merge(y    = pData(eset),
          by.x = "rowname",
          by.y = "uid")

  if(returnObject == "probSubjects"){
    colors <- RColorBrewer::brewer.pal(12, "Paired")
    set.seed(10)
    colors <- sample(colors, length(colors))
    colors <- grDevices::colorRampPalette(colors)(length(unique(plotDF$participant_id)))

    fullPlot <- ggplot2::ggplot(data = plotDF,
                                mapping = ggplot2::aes(x = y_chrom_present, y = chry)) +
      ggplot2::geom_boxplot(outlier.color = "transparent", fill = "grey") +
      ggplot2::geom_jitter(height = 0, width = 0.25, mapping = ggplot2::aes(color = participant_id)) +
      ggplot2::labs(y = "Average probe intensities (chrY)") +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::facet_wrap(facets = ~study_accession+matrix, scale = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.pos = "none",
                     strip.text = ggplot2::element_text(size = 5))
    return(fullPlot)
  }

  # find outlier (1.5 x IGR from Q1 and Q3)
  outlierDF <- plotDF %>%
    dplyr::group_by(study_accession, matrix, y_chrom_present) %>%
    dplyr::mutate(up = chry >= stats::quantile(chry, probs = 0.75) + 1.5 * stats::IQR(chry),
           dn = chry <= stats::quantile(chry, probs = 0.25) - 1.5 * stats::IQR(chry))

  quartileDF <- plotDF %>%
    dplyr::group_by(study_accession, matrix, y_chrom_present) %>%
    dplyr::summarize(q1 = stats::quantile(chry, probs = 0.25) - 1.5 * stats::IQR(chry),
              q3 = stats::quantile(chry, probs = 0.75) + 1.5 * stats::IQR(chry)) %>%
    dplyr::mutate(y_chrom_present = !y_chrom_present)

  # flag outlier samples (possible swap)
  flagDF <- merge(x  = outlierDF,
                  y  = quartileDF,
                  by = c("study_accession", "matrix", "y_chrom_present")) %>%
    dplyr::mutate(flag = FALSE,
           flag = ifelse(test = !y_chrom_present & up & chry >= q1,
                         yes  = TRUE,
                         no   = flag),
           flag = ifelse(test = y_chrom_present & dn & chry <= q3,
                         yes  = TRUE,
                         no   = flag))

  plotDF <- merge(x     = plotDF,
                  y     = dplyr::select(flagDF, rowname, flag),
                  by    = "rowname",
                  all.x = TRUE) %>%
    # if y_chrom_present not specified, set flag as false
    dplyr::mutate(flag = ifelse(test = is.na(flag),
                         yes  = FALSE,
                         no   = flag))



  if(returnObject == "allMatricesPlot"){
    fullPlot <- ggplot2::ggplot(data = plotDF,
                       mapping = ggplot2::aes(x = y_chrom_present, y = chry)) +
      ggplot2::geom_boxplot(outlier.color = "transparent", fill = "grey") +
      ggplot2::geom_jitter(height = 0, width = 0.25, mapping = ggplot2::aes(color = flag)) +
      ggplot2::labs(y = "Average probe intensities (chrY)") +
      ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      ggplot2::facet_wrap(facets = ~study_accession+matrix, scale = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
            legend.pos = "none",
            strip.text = ggplot2::element_text(size = 5))
    return(fullPlot)

  }else if(returnObject == "probSamplesDT"){
    return(plotDF[ plotDF$flag, ])

  }else if(returnObject == "probSamplesPlot"){
    problemMatrix <- unique(plotDF$matrix[ plotDF$flag])
    plotTemp <- dplyr::filter(plotDF, matrix %in% problemMatrix)

    prbSmplsPlot <- ggplot2::ggplot(data = plotTemp,
                           mapping = ggplot2::aes(x = y_chrom_present, y = chry)) +
      ggplot2::geom_boxplot(outlier.color = "transparent", fill = "grey") +
      ggbeeswarm::geom_beeswarm(mapping = ggplot2::aes(color = flag), cex = 2.5, size = 0.8) +
      ggplot2::labs(y = "Average probe intensities (chrY)", x = "Sex") +
      ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      ggplot2::facet_wrap(facets = ~study_accession+matrix, scale = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
            legend.pos = "none",
            strip.text = ggplot2::element_text(size = 6))
    return(prbSmplsPlot)
  }else{
    stop("returnObject parameter value not recognized")

  }
}

#' Generate gender box plots (prior to QC)
#'
#' @param eset expressionSet
#' @export
#'
qualityControl.genderByMatrix <- function(eset){
  yChromGenes <- yChromGenes[ yChromGenes %in% featureNames(eset) ]

  plotDF <- colMeans(exprs(eset)[yChromGenes, ], na.rm = TRUE) %>%
    data.frame(chry = .) %>%
    tibble::rownames_to_column() %>%
    merge(y    = pData(eset),
          by.x = "rowname",
          by.y = "uid")

  outlierDF <- plotDF %>%
    dplyr::group_by(study_accession, matrix, gender) %>%
    dplyr::mutate(up = chry >= stats::quantile(chry, probs = 0.75) + 1.5 * stats::IQR(chry),
           dn = chry <= stats::quantile(chry, probs = 0.25) - 1.5 * stats::IQR(chry))

  quartileDF <- plotDF %>%
    dplyr::group_by(study_accession, matrix, gender) %>%
    dplyr::summarize(q1 = stats::quantile(chry, probs = 0.25) - 1.5 * stats::IQR(chry),
                     q3 = stats::quantile(chry, probs = 0.75) + 1.5 * stats::IQR(chry)) %>%
    dplyr::mutate(gender = c("Female" = "Male", "Male" = "Female")[gender])

  # flag outlier samples (possible swap)
  flagDF <- merge(x  = outlierDF,
                  y  = quartileDF,
                  by = c("study_accession", "matrix", "gender")) %>%
    dplyr::mutate(flag = FALSE,
                  flag = ifelse(test = gender %in% "Female" & up & chry >= q1,
                                yes  = TRUE,
                                no   = flag),
                  flag = ifelse(test = gender %in% "Male" & dn & chry <= q3,
                                yes  = TRUE,
                                no   = flag))

  plotDF <- merge(x     = plotDF,
                  y     = dplyr::select(flagDF, rowname, flag),
                  by    = "rowname",
                  all.x = TRUE) %>%
    # if gender not specified, set flag as false
    dplyr::mutate(flag = ifelse(test = is.na(flag),
                                yes  = FALSE,
                                no   = flag))

  fullPlot <- ggplot2::ggplot(data = plotDF,
                              mapping = ggplot2::aes(x = gender, y = chry)) +
    ggplot2::geom_boxplot(outlier.color = "transparent", fill = "grey") +
    ggplot2::geom_jitter(height = 0, width = 0.25, mapping = ggplot2::aes(color = flag)) +
    ggplot2::labs(y = "Average probe intensities (chrY)") +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    ggplot2::facet_wrap(facets = ~study_accession+matrix, scale = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.pos = "none",
                   strip.text = ggplot2::element_text(size = 5))
}

#' Generate ychrom box plots after QC
#'
#' @param eset expressionSet
#' @export
#'
qualityControl.failedYchromQC <- function(eset){
  yChromGenes <- yChromGenes[ yChromGenes %in% featureNames(eset) ]

  plotDF <- colMeans(exprs(eset)[yChromGenes, ], na.rm = TRUE) %>%
    data.frame(chry = .) %>%
    tibble::rownames_to_column()

  plotDF <- merge(plotDF, pData(eset), by.x = "rowname", by.y = "uid")

  outlierDF <- plotDF %>%
    dplyr::group_by(study_accession, matrix, y_chrom_present) %>%
    dplyr::mutate(up = chry >= stats::quantile(chry, probs = 0.75) + 1.5 * stats::IQR(chry),
                  dn = chry <= stats::quantile(chry, probs = 0.25) - 1.5 * stats::IQR(chry))

  quartileDF <- plotDF %>%
    dplyr::group_by(study_accession, matrix, y_chrom_present) %>%
    dplyr::summarize(q1 = stats::quantile(chry, probs = 0.25) - 1.5 * stats::IQR(chry),
                     q3 = stats::quantile(chry, probs = 0.75) + 1.5 * stats::IQR(chry)) %>%
    dplyr::mutate(y_chrom_present = !y_chrom_present)

  # flag outlier samples (possible swap)
  flagDF <- merge(x  = outlierDF,
                  y  = quartileDF,
                  by = c("study_accession", "matrix", "y_chrom_present"))

  plotDF <- merge(x     = plotDF,
                  y     = dplyr::select(flagDF, rowname),
                  by    = "rowname",
                  all.x = TRUE)

  fullPlot <- ggplot2::ggplot(data = plotDF,
                            mapping = ggplot2::aes(x = y_chrom_present, y = chry)) +
    ggplot2::geom_boxplot(outlier.color = "transparent", fill = "grey") +
    ggplot2::geom_jitter(height = 0, width = 0.25, mapping = ggplot2::aes(color = failedYchromQC)) +
    ggplot2::labs(y = "Average probe intensities (chrY)") +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    ggplot2::facet_wrap(facets = ~study_accession+matrix, scale = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                 legend.pos = "none",
                 strip.text = ggplot2::element_text(size = 5))
}
#' Generate table of studies with count of participants from vector of pids
#'
#' @param participantIds participant IDs
#' @export
#'
qualityControl.createSubjectsByStudyTable <- function(participantIds){
  uniquePids <- unique(participantIds)
  studies <- gsub("SUB\\d{6}\\.", "SDY", uniquePids)
  studiesTbl <- table(studies)
}
