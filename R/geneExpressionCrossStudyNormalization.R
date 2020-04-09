#' Cross-Study normalize an expressionSet object based on a single
#' vendor's subset of expression data to define the target distribution
#'
#' @param eset expressionSet
#' @param vendorToUse Microarray or RNAseq vendor to use for generating target distribution
#' @param studiesToExclude studies to exclude when generating the target distribution
#' @export
#'
crossStudyNormalize <- function(eset, vendorToUse, studiesToExclude){
  # Find samples having at least 60% of genes with non-NA value
  geneSymbolMeans <- colMeans(is.na(exprs(eset)))
  samplesToUse <- names(geneSymbolMeans)[ geneSymbolMeans > 0.6]

  # Find features that are present in samples with good coverage of gene expression
  tmp <- exprs(eset[, samplesToUse])
  featuresToUse <- rownames(tmp)[ complete.cases(tmp) ]
  rm(tmp)

  # Use only Affymetrix samples to generate target distribution as this ensures
  # a reasonable range of values for all studies to use as it removes cross-platform
  # technical variation in range of expression values which causes the distribution
  # to be too wide.  SDY1293 also generated a flattened distribution and was excluded.
  exprs.not.norm <- exprs(eset)
  excludedSmpls <- which(eset$featureSetVendor != vendorToUse | eset$study_accession %in% studiesToExclude)
  targetExprs <- exprs.not.norm[ featuresToUse, -excludedSmpls]

  # perform normalization on selected features (12127) and samples (1654)
  normTargetExprs <- preprocessCore::normalize.quantiles.robust(targetExprs)

  # Create target distribution based on normalized version of selected features and samples
  target.dist <- preprocessCore::normalize.quantiles.determine.target(normTargetExprs)

  # Clean up as we go
  rm(normTargetExprs, targetExprs)

  # Normalize all features and samples based on subset target distribution
  # Even though selected features are only 12127 and the number of features per sample may be more or less.
  # It is unclear how this problem is handled.
  normAllExprs <- preprocessCore::normalize.quantiles.use.target(exprs.not.norm, target = target.dist)

  # replace expression values in expressionSet and create copy with clear name
  dimnames(normAllExprs) <- dimnames(exprs.not.norm)
  exprs(eset) <- normAllExprs

  # Clean up as we go
  rm(normAllExprs)
  gc()

  return(eset)
}

#' Create new featureSetName field that coalesces similar platforms for analysis
#'
#' @param eset expressionSet
#' @export
#'
addCoalescedFeatureSetName <- function(eset){
  eset$featureSetName2 <- eset$featureSetName
  eset$featureSetName2[ grepl('ustom', eset$featureSetName) ] <- 'RNA-seq'
  eset$featureSetName2[ grepl('HT-12', eset$featureSetName) ] <- 'HumanHT-12_2018'
  return(eset)
}

#' Adjust baseline expression data to account for study_accession, featureSetName, and cell_type
#'
#' @param eset expressionSet
#' @export
#'
batchCorrectBaselineData <- function(eset){
  day0 <- eset[, eset$time_post_last_vax == 0]

  # Create model matrix for selected variables in baseline data eset
  model.vars <- c('cell_type','gender_imputed','featureSetName2','study_accession2')
  mm <- model.matrix(as.formula(paste0('~',paste0(model.vars, collapse='+'))),
                     data = pData(day0))

  # Remove model matrix variables that are not estimable
  # unique FAS: study2SDY1291 study2SDY1289 study2SDY1370 study2SDY80
  # ???: study2SDY1294 study2SDY61 study2SDY63 study2SDY67
  notEstimable <- nonEstimable(mm)
  mm.est <- mm[, !colnames(mm) %in% notEstimable]

  # Create a linear model based on the model matrix with all estimable vars for day 0
  fit <- lmFit(object = exprs(day0), design = mm.est)

  # Select only study, featureAnnotationSet, and cell type as variables to use
  # (do not use gender)
  coefs2adjust <- grep('study|featureSetName|cell_type', colnames(mm.est), value = TRUE)

  # create another model matrix with all estimable vars for all timepoints
  mm.all <- model.matrix(as.formula(paste0('~',paste0(model.vars, collapse='+'))),
                         data = pData(eset))
  mm.all <- mm.all[,(colnames(mm.all) %in% colnames(mm.est))]

  # subset model matrix of all timepoints and estimable variables by the coefficients to adjust
  mm.subset <- t(mm.all[, coefs2adjust])

  # subset the coefficient values of the linear fit output by coefficients to adjust as well
  fit.vals <- fit$coefficients[, coefs2adjust]

  # adjust the expression values at baseline for each sample and feature by an amount proportional to the sum of the coefficient values multiplied by the presence or absence of the coefficient variable
  exprs(eset) <- exprs(eset) - fit.vals %*% mm.subset

  return(eset)
}

#' Generate sample MDS plots for QC
#'
#' @param eset expressionSet
#' @export
#'
generateSampleMDSPlot <- function(eset, numberOfSamples){
  # Based on code by Aris - 03/10/2019
  exp <- exprs(eset)
  pheno <- pData(eset)

  # Get data for time point = 0; remove NA genes
  dset <- exp[, rownames(pheno)[pheno$time_post_last_vax == 0]]
  keep <- apply(dset, 1, function(x){return(sum(is.na(x)) == 0)})
  dset <- dset[keep, ]

  # To make plotting more clear, use hardcoded number of samples from each study
  groups <- pheno[colnames(dset), ]$study_accession
  names(groups) <- colnames(dset)
  ind <- unlist(lapply(unique(groups), function(x){
    smpls <- which(groups == x)
    toReturn <- ifelse(length(smpls) > numberOfSamples, numberOfSamples, length(smpls))
    return(smpls[1:toReturn])
  }))
  dset <- dset[, ind]
  groups <- groups[colnames(dset)]
  colors <- RColorBrewer::brewer.pal(10, "Spectral")
  colors <- colorRampPalette(colors)(length(unique(groups)))
  tmp <- plotMDS(dset, col = colors[factor(groups)], labels = groups)
}

#' Remove incomplete rows
#'
#' @param eset expressionSet
#' @export
#'
removeAllNArows <- function(eset){
  em <- Biobase::exprs(eset)
  allNARows <- apply(em, 1, function(x){ all(is.na(x))})
  eset <- eset[ !allNARows, ]
}
