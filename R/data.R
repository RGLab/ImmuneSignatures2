
NULL

#' Vaccine Mapping Data.
#'
#' Mapping of vaccine meta-data to gene expression matrices created by
#' Thomas Hagan and Evan Henrich
#'
#' @format A data frame with seven variables: \code{matrix}, \code{vaccine},
#'   \code{vaccine_type}, \code{pathogen}, \code{adjuvant}, \code{study}
#'   and \code{arm_accesion}
"vaccines"

#' SDY1370 ELISA data
#'
#' SDY1370 author provided ELISA results formatted to work with ImmuneSpace-derived data
#'
#' @format A data frame that mimics ImmuneSpace `con$getDataset("elisa", original_view=TRUE)` results
"sdy1370_elisa"

#' SDY1325 meta-data
#
#' Derived from study-associated paper's figure 1.
#'
#' @format data.table
"sdy1325_metadata"

#' Y-Chromosome Associated Genes
#'
#' Y-Chromosome genes from the ensembl database using the `biomaRt` library
#'
#' @format string vector containing all gene names
"yChromGenes"

#' Expected Gene Expression Meta-Data columns
#'
#' After pre-processing, a vector of column names expected to be found in geMetaData
#'
#' @format string vector containing all column names
"expectedGeMetaDataColumns"

#' Age groups for determining age cohorts
#'
#' Data created from inputSamples query in ImmuneSpace. Age_bracket variable
#' is the ceiling of that decile.
#'
#' @format A data frame with three variables: \code{count}, \code{study},
#'   and \code{age_bracket}
"ageGroups"
