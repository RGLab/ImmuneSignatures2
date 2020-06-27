#' Cross-Study normalize an expressionSet object based on a single
#' vendor's subset of expression data to define the target distribution
#'
#' @param eset expressionSet
#' @param targetDistributionVendor Microarray or RNAseq vendor to use for generating target distribution
#' @param targetDistributionExcludedStudies studies to exclude when generating the target distribution
#' @export
#'
crossStudyNormalize <- function(eset, targetDistributionVendor, targetDistributionExcludedStudies){
  # Find samples having at least 50% of genes with non-NA value ... necessary?
  geneSymbolMeans <- colMeans(is.na(exprs(eset)))
  samplesToUse <- names(geneSymbolMeans)[ geneSymbolMeans > 0.5 ]
  if(length(samplesToUse) == 0){
    warning("no samples found with sufficient coverage ... using all")
  }

  # Find features that are present in samples with good coverage of gene expression
  tmp <- exprs(eset[, samplesToUse])
  featuresToUse <- rownames(tmp)[ complete.cases(tmp) ]
  rm(tmp)

  # Use only Affymetrix samples to generate target distribution as this ensures
  # a reasonable range of values for all studies to use as it removes cross-platform
  # technical variation in range of expression values which causes the distribution
  # to be too wide.  SDY1293 also generated a flattened distribution and was excluded.
  exprs.not.norm <- exprs(eset)
  excludedSmpls <- which(eset$featureSetVendor != targetDistributionVendor |
                         eset$study_accession %in% targetDistributionExcludedStudies)
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

  # Ensure formatting of returned eSet
  pd <- pData(eset)
  rownames(pd) <- pd$uid
  normAllExprs <- normAllExprs[ , order(match(colnames(normAllExprs), rownames(pd))) ]
  normEset <- new("ExpressionSet",
                  exprs = as.matrix(normAllExprs, rownames = rownames(normAllExprs)),
                  phenoData = new('AnnotatedDataFrame',
                                  pd)
  )

  # Clean up as we go
  rm(normAllExprs)
  gc()

  return(normEset)
}

#' Adjust baseline expression data to account for study_accession, featureSetName, and cell_type
#'
#' @param eset expressionSet
#' @export
#'
batchCorrectBaselineData <- function(eset){
  baseline <- eset[, eset$time_post_last_vax <= 0]

  # Create model matrix for selected variables in baseline data eset
  model.vars <- c('gender_imputed','cell_type','featureSetName2','geBatchName')
  mm <- model.matrix(as.formula(paste0('~',paste0(model.vars, collapse='+'))),
                     data = pData(baseline))

  # Remove model matrix variables that are not estimable
  # unique FAS: study2SDY1291 study2SDY1289 study2SDY1370 study2SDY80
  # ???: study2SDY1294 study2SDY61 study2SDY63 study2SDY67
  notEstimable <- nonEstimable(mm)
  mm.est <- mm[, !colnames(mm) %in% notEstimable]

  # Create a linear model based on the model matrix with all estimable vars for day 0
  fit <- lmFit(object = exprs(baseline), design = mm.est)

  # Select only study, featureAnnotationSet, and cell type as variables to use
  # (do not use gender)
  coefs2adjust <- grep('Batch|featureSetName|cell_type', colnames(mm.est), value = TRUE)

  # create another model matrix with all estimable vars for all timepoints
  mm.all <- model.matrix(as.formula(paste0('~', paste0(model.vars, collapse='+'))),
                         data = pData(eset))
  mm.all <- mm.all[,(colnames(mm.all) %in% colnames(mm.est))]

  # subset model matrix of all timepoints and estimable variables by the coefficients to adjust
  mm.subset <- t(mm.all[, coefs2adjust])

  # subset the coefficient values of the linear fit output by coefficients to adjust as well
  fit.vals <- fit$coefficients[, coefs2adjust]

  # adjust the expression values at baseline for each sample and feature by an amount proportional to the sum of the coefficient values multiplied by the presence or absence of the coefficient variable
  correctionValues <- fit.vals %*% mm.subset
  exprs(eset) <- exprs(eset) - correctionValues

  return(eset)
}


