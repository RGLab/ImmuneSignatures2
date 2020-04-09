
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

#' Y-Chromosome Associated Genes
#'
#' Y-Chromosome genes from the ensembl database using the `biomaRt` library
#'
#' @format string vector containing all gene names
"yChromGenes"

#' Age groups for determining age cohorts
#'
#' Data created from inputSamples query in ImmuneSpace. Age_bracket variable
#' is the ceiling of that decile.
#'
#' @format A data frame with three variables: \code{count}, \code{study},
#'   and \code{age_bracket}
"ageGroups"
