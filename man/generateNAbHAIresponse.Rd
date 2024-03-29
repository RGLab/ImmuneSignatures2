% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immuneResponseCallGeneration.R
\name{generateNAbHAIresponse}
\alias{generateNAbHAIresponse}
\title{Generate immune response calls for HAI or NAb using pipeline originally
developed by Daniel Chawla at Yale University}
\usage{
generateNAbHAIresponse(assay, df, postVaxDayRange, discretizationValues)
}
\arguments{
\item{assay}{assay name}

\item{df}{immdata df}

\item{postVaxDayRange}{Allowable timepoints for post-vaccine values}

\item{discretizationValues}{cut points to use for discretizing response call groups}
}
\description{
Generate immune response calls for HAI or NAb using pipeline originally
developed by Daniel Chawla at Yale University
}
