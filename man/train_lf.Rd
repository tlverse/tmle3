% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{train_lf}
\alias{train_lf}
\title{Manually Train Likelihood Factor
The internal training process for likelihood factors is somewhat obtuse, so this function
does the steps to manually train one, which is helpful if you want to use a likelihood factor
independently of a likelihood object}
\usage{
train_lf(lf, tmle_task)
}
\arguments{
\item{lf}{the likelihood factor to train}

\item{tmle_task}{the task to use for training}
}
\description{
Manually Train Likelihood Factor
The internal training process for likelihood factors is somewhat obtuse, so this function
does the steps to manually train one, which is helpful if you want to use a likelihood factor
independently of a likelihood object
}
