% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmle3_Spec_stratified.R
\name{tmle_stratified}
\alias{tmle_stratified}
\title{Stratified version of TML estimator from other Spec classes}
\usage{
tmle_stratified(base_spec, strata_variable)
}
\arguments{
\item{base_spec}{An underlying spec to stratify.}

\item{strata_variable}{The variable(s) to use for stratification.}
}
\description{
O=(W,A,Y)
W=Covariates
A=Treatment (binary or categorical)
Y=Outcome (binary or bounded continuous)
}
