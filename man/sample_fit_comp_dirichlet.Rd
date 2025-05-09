% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_comp_dirichlet.R
\name{sample_fit_comp_dirichlet}
\alias{sample_fit_comp_dirichlet}
\title{Sample data and fit a Dirichlet regression model}
\usage{
sample_fit_comp_dirichlet(
  data,
  n,
  method,
  bias_var = NULL,
  bias_thresh = NULL,
  clh_var = NULL,
  clh_iter = NULL,
  strata_var = NULL,
  formula
)
}
\arguments{
\item{data}{A data frame or a \code{PackedSpatRaster} object containing the input data.}

\item{n}{Integer. The number of samples to draw.}

\item{method}{Character. The sampling method to use (e.g., "random", "stratified").}

\item{bias_var}{Character. The name of the variable to use for biased sampling (optional).}

\item{bias_thresh}{Numeric. A threshold for biased sampling (optional).}

\item{clh_var}{Character. The name of the variable for conditional Latin hypercube sampling (optional).}

\item{clh_iter}{Integer. The number of iterations for conditional Latin hypercube sampling (optional).}

\item{strata_var}{Character. The name of the variable to use for stratified sampling (optional).}

\item{formula}{A formula specifying the regression model to be fitted.}
}
\value{
A modified version of the input \code{data} object, including model coefficients
and predicted substrate composition probabilities.
}
\description{
This function samples the input data and fits a Dirichlet regression model using
the \code{fit_comp_dirichlet()} function.
}
