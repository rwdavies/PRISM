% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{prepRepeat}
\alias{prepRepeat}
\title{Prep repeat masker file}
\usage{
prepRepeat(
  rmsk,
  simpleRepeats,
  chrList,
  verbose = 1,
  nCores = 1,
  useChr = TRUE
)
}
\arguments{
\item{rmsk}{Path to repeat masker file}

\item{simpleRepeats}{Path to simple repeats file}

\item{chrList}{vector of chromosome names to prep}

\item{verbose}{whether to print messages. verbose 1 = Yes, 0 = No}

\item{nCores}{number of cores to us}

\item{useChr}{Whether to add chr to repeat masker genoName list}
}
\description{
Prep repeat masker file
}
