#' Cross-Study normalize an expressionSet object based on a single
#' vendor's subset of expression data to define the target distribution
#'
#' @param eset expressionSet
#' @param targetDistributionVendor Microarray or RNAseq vendor to use for generating target distribution
#' @param targetDistributionExcludedStudies studies to exclude when generating the target distribution
#' @export
#'
crossStudyNormalize <- function(eset, targetDistributionVendor, targetDistributionExcludedStudies){
  eset.target <- eset[ , eset$featureSetVendor == targetDistributionVendor &
                         !(eset$study_accession %in% targetDistributionExcludedStudies) ]

  # Ensure we only use features with complete.cases, as normalize.quantiles.robust
  # doesn't handle NAs well - shrinks quantiles
  eset.target <- eset.target[ stats::complete.cases(exprs(eset.target)), ]

  # perform normalization on selected features (14827) and samples (1973)
  normTargetExprs <- preprocessCore::normalize.quantiles.robust(exprs(eset.target), use.log2 = TRUE)

  # Create target distribution based on normalized version of selected features and samples
  target.dist <- preprocessCore::normalize.quantiles.determine.target(normTargetExprs)

  # Normalize all features and samples based on subset target distribution
  normAllExprs <- preprocessCore::normalize.quantiles.use.target(exprs(eset), target = target.dist)

  # replace expression values in expressionSet and create copy with clear name
  dimnames(normAllExprs) <- dimnames(exprs(eset))

  # Ensure formatting of returned eSet
  pd <- pData(eset)
  rownames(pd) <- pd$uid
  normAllExprs <- normAllExprs[ , order(match(colnames(normAllExprs), rownames(pd))) ]
  normEset <- methods::new("ExpressionSet",
                  exprs = as.matrix(normAllExprs, rownames = rownames(normAllExprs)),
                  phenoData = methods::new('AnnotatedDataFrame',
                                  pd)
  )

  # Clean up as we go
  rm(normAllExprs, normTargetExprs, eset.target)
  gc()

  return(normEset)
}

#' Adjust expression data using correction coefficients based on baseline data
#' to account for study_accession, featureSetName, and cell_type
#' @param eset expressionSet
#' @param batch.vars variables to correct for
#' @param covariates additional variables to include in the model
#' @export
#'
batchCorrect <- function(eset, batch.vars, covariates){
  eset.baseline <- eset[, eset$time_post_last_vax <= 0 ]

  model.vars <- c(batch.vars, covariates)
  # samples that failed Y-chromosome QC are not used for generating coefficients,
  # but are left in the dataset for analysts to use / rm
  eset.baseline.passedYchromQC <- eset.baseline[ , !eset.baseline$failedYchromQC ]

  model.formula <- getModelFormula(model.vars)

  fit.vals <- getLmFitValues(model.formula,
                             eset.baseline.passedYchromQC,
                             coefs2adjust = batch.vars)

  # create another model matrix with all estimable vars for all timepoints
  mm.all <- stats::model.matrix(model.formula, data = pData(eset))
  mm.all.subset <- mm.all[, (colnames(mm.all) %in% colnames(fit.vals))]

  # adjust the expression values for each sample and feature
  # by an amount proportional to the sum of the coefficient values multiplied
  # by the presence or absence of the coefficient variable
  correctionValues <- fit.vals %*% t(mm.all.subset)
  tmpExprs <- exprs(eset) - correctionValues

  eset.corr <- methods::new("ExpressionSet",
                   exprs = as.matrix(tmpExprs, rownames = rownames(tmpExprs)),
                   phenoData = methods::new('AnnotatedDataFrame', pData(eset))
  )
  return(eset.corr)
}

#' Adjust baseline expression data to account for study_accession, featureSetName, and cell_type
#'
#' @param modelEset expressionSet used to create linear model
#' @param targetEset expressionSet that is corrected given model fit values
#' @param batch.vars variables to correct for
#' @param covariates additional variables to include in the model
#' @export
#'
batchCorrect.importedModel <- function(modelEset, targetEset, batch.vars, covariates){
  model.vars <- c(batch.vars, covariates)
  esets <- c(model = modelEset, target = targetEset)
  factorsToRelevel <- c(model.vars)
  # Make sure "young" comes first
  esets <- relevelEsets(esets, factorsToRelevel)

  modelEset <- esets$model
  targetEset <- esets$target

  # Ensure baseline and passing Y-chromosome QC only for determining coefficients
  modelEset.baseline <- modelEset[, modelEset$time_post_last_vax <= 0 ]
  modelEset.baseline.passedYchromQC <- modelEset.baseline[, !modelEset.baseline$failedYchromQC ]

  model.formula <- getModelFormula(model.vars)
  fit.vals <- getLmFitValues(model.formula,
                             modelEset.baseline.passedYchromQC,
                             coefs2adjust =  batch.vars)

  exprs.target <- exprs(targetEset)
  geneIntersect <- intersect(rownames(exprs.target),
                             rownames(exprs(modelEset.baseline.passedYchromQC)))
  exprs.target.subset <- exprs.target[ geneIntersect, ]

  # Find common components of the model matrix for target and the modelEset.baseline
  mm.target <- stats::model.matrix(model.formula, data = pData(targetEset))
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


  # ensure phenotypic data is ordered correctly as this is not checked
  # when inserted directly - i.e. pData(x) <- pd
  pd <- pData(targetEset)
  pd <- pd[ order(match(pd$uid, colnames(exprs.target.corr))), ]

  # Create new expressionSet object as validation fails for insertion into target
  # due to internal checks related to expression matrix
  eset.target.corr <- methods::new("ExpressionSet",
                       exprs = as.matrix(exprs.target.corr, rownames = rownames(exprs.target.corr)),
                       phenoData = methods::new('AnnotatedDataFrame', pd))
}

#' Get linear model fit values for each gene
#'
#' @param eset.baseline baseline expressionSet
#' @param model.formula model formula
#' @param coefs2adjust Coeffieicts to include in the adjusted model
#' @export
#'
getLmFitValues <- function(model.formula, eset.baseline, coefs2adjust){
  mm <- stats::model.matrix(model.formula, data = pData(eset.baseline))

  # Remove model matrix variables that are not estimable
  notEstimable <- limma::nonEstimable(mm)
  mm.est <- mm[, !colnames(mm) %in% notEstimable]

  # Create a linear model based on the model matrix with all estimable vars for day 0
  fit <- limma::lmFit(object = exprs(eset.baseline), design = mm.est)

  # Select only study, featureAnnotationSet, and cell type as variables to use
  # (do not use gender)
  coefs2adjust <- grep(paste0(coefs2adjust, collapse = "|"),
                       colnames(mm.est), value = TRUE)

  # subset the coefficient values of the linear fit output by coefficients to adjust as well
  fit.vals <- fit$coefficients[, coefs2adjust]
}

#' Get linear model formula
#'
#' @param model.vars vector of variables
#' @export
#'
getModelFormula <- function(model.vars){
  model.formula <- stats::as.formula(paste0('~', paste0(model.vars, collapse='+')))
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
    esets$model[[fac]] <- stats::relevel(as.factor(esets$model[[fac]]),
                                           ref = ref.batch)
    esets$target[[fac]] <- stats::relevel(as.factor(esets$target[[fac]]),
                                         ref = ref.batch)
  }
  return(esets)
}
