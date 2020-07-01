#' Cross-Study normalize an expressionSet object based on a single
#' vendor's subset of expression data to define the target distribution
#'
#' @param eset expressionSet
#' @param targetDistributionVendor Microarray or RNAseq vendor to use for generating target distribution
#' @param targetDistributionExcludedStudies studies to exclude when generating the target distribution
#' @export
#'
crossStudyNormalize <- function(eset, targetDistributionVendor, targetDistributionExcludedStudies){
  # Find samples having at least 60% of genes with non-NA value ... necessary?
  geneSymbolMeans <- colMeans(!is.na(exprs(eset)))
  samplesToUse <- names(geneSymbolMeans)[ geneSymbolMeans > 0.6 ]
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
batchCorrect <- function(eset){
  eset.baseline <- eset[, eset$time_post_last_vax <= 0 ]

  model.formula <- getModelFormula()

  fit.vals <- getLmFitValues(model.formula, eset.baseline)

  # create another model matrix with all estimable vars for all timepoints
  mm.all <- model.matrix(model.formula, data = pData(eset))
  mm.all.subset <- mm.all[, (colnames(mm.all) %in% colnames(fit.vals))]

  # adjust the expression values at baseline for each sample and feature
  # by an amount proportional to the sum of the coefficient values multiplied
  # by the presence or absence of the coefficient variable
  correctionValues <- fit.vals %*% t(mm.all.subset)
  tmpExprs <- exprs(eset) - correctionValues

  eset.corr <- new("ExpressionSet",
                   exprs = as.matrix(tmpExprs, rownames = rownames(tmpExprs)),
                   phenoData = new('AnnotatedDataFrame', pData(eset))
  )

}

#' Adjust baseline expression data to account for study_accession, featureSetName, and cell_type
#'
#' @param modelEset expressionSet used to create linear model
#' @param targetEset expressionSet that is corrected given model fit values
#' @export
#'
batchCorrect.importedModel <- function(modelEset, targetEset){
  esets <- c(model = modelEset, target = targetEset)
  factorsToRelevel <- c("cell_type",
                        "featureSetVendor",
                        "geBatchName",
                        "featureSetName")
  esets <- relevelEsets(esets, factorsToRelevel)

  modelEset <- esets$model
  targetEset <- esets$target

  modelEset.baseline <- modelEset[, modelEset$time_post_last_vax <= 0]

  model.formula <- getModelFormula()
  fit.vals <- getLmFitValues(model.formula, modelEset.baseline)

  exprs.target <- exprs(targetEset)
  geneIntersect <- intersect(rownames(exprs.target),
                             rownames(exprs(modelEset.baseline)))
  exprs.target.subset <- exprs.target[ geneIntersect, ]

  # Find common components of the model matrix for target and the modelEset.baseline
  mm.target <- model.matrix(model.formula, data = pData(targetEset))
  factorIntersect <- intersect(colnames(mm.target), colnames(fit.vals))
  mm.target.subset <- mm.target[ , factorIntersect ]

  # Subset fit values to match dimensions of
  fit.vals.subset <- fit.vals[ geneIntersect, factorIntersect ]

  # Create correction values by matrix multiplication of the gene * coef fit matrix
  mm.target.subset.t <- t(mm.target.subset)

  # Double-check the column order for yfit-vals and mm.target.subset
  mm.target.subset.t <- mm.target.subset.t[ order(match(rownames(mm.target.subset.t),
                                                  colnames(fit.vals.subset))), ]

  # Ensure that correction matrix has same ordering as the expression matrix to be corrected
  correctionValues <- fit.vals.subset %*% mm.target.subset.t
  correctionValues <- correctionValues[ order(match(rownames(correctionValues),
                                                    rownames(exprs.target.subset))), ]

  # create new expression matrix
  exprs.target.corr <- exprs.target.subset - correctionValues

  # pd
  pd <- pData(targetEset)
  pd <- pd[ order(match(pd$uid, colnames(exprs.target.corr))), ]

  # Create new expressionSet object as validation fails for insertion into target
  eset.target.corr <- new("ExpressionSet",
                       exprs = as.matrix(exprs.target.corr, rownames = rownames(exprs.target.corr)),
                       phenoData = new('AnnotatedDataFrame', pd))
}

#' Get linear model fit values for each gene
#'
#' @param eset.baseline baseline expressionSet
#' @param model.formula model formula
#' @export
#'
getLmFitValues <- function(model.formula, eset.baseline){
  mm <- model.matrix(model.formula, data = pData(eset.baseline))

  # Remove model matrix variables that are not estimable
  notEstimable <- nonEstimable(mm)
  mm.est <- mm[, !colnames(mm) %in% notEstimable]

  # Create a linear model based on the model matrix with all estimable vars for day 0
  fit <- lmFit(object = exprs(eset.baseline), design = mm.est)

  # Select only study, featureAnnotationSet, and cell type as variables to use
  # (do not use gender)
  coefs2adjust <- grep('Batch|featureSetName|featureSetVendor|cell_type',
                       colnames(mm.est), value = TRUE)

  # subset the coefficient values of the linear fit output by coefficients to adjust as well
  fit.vals <- fit$coefficients[, coefs2adjust]
}

#' Get linear model fit values for each gene
#'
#' @export
#'
getModelFormula <- function(){
  model.vars <- c('gender_imputed','cell_type','featureSetName','geBatchName','featureSetVendor')
  model.formula <- as.formula(paste0('~', paste0(model.vars, collapse='+')))
}

#' Get linear model fit values for each gene
#'
#' @param esets expressionSets
#' @param factorsToRelevel factors to relevel
#' @export
#'
relevelEsets <- function(esets, factorsToRelevel){
  for(fac in factorsToRelevel){
    ref.batch <- intersect(unique(esets$model[[fac]]),
                           unique(esets$target[[fac]]))[[1]]
    esets$model[[fac]] <- relevel(as.factor(esets$model[[fac]]),
                                           ref = ref.batch)
    esets$target[[fac]] <- relevel(as.factor(esets$target[[fac]]),
                                         ref = ref.batch)
  }
  return(esets)
}
