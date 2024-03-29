% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_vaxdt.R
\name{derive_vars_vaxdt}
\alias{derive_vars_vaxdt}
\title{Add Vaccination Date Variables to the Output Dataset}
\usage{
derive_vars_vaxdt(dataset, dataset_adsl, by_vars, order)
}
\arguments{
\item{dataset}{Input dataset

The variables specified by the \code{by_vars} argument are expected.}

\item{dataset_adsl}{Input adsl dataset

The vaccination date variables created will be merged with this adsl dataset.}

\item{by_vars}{Grouping variables.

The variables to be grouped to filter the first observation within each
by group.}

\item{order}{Sorting variables.

The variables order to be specified either in ascending or descending order.
By default ascending order will be applicable.}
}
\value{
The adsl dataset with vaccination date variables added to it.
}
\description{
Creates vaccination date variables from \code{EX} domain. A date variable will be
created for each vaccination taking values from the variable \code{EXSTDTC}.
}
\details{
If there are multiple vaccinations for a visit per subject,warning will be
provided and only first observation will be filtered based on the variable
order specified on the \code{order} argument. In this case, user need to select
the \code{by_vars} appropriately.

The number of variables created will be based on the number of vaccinations
per subject per visit.
}
\examples{
library(tibble)
library(admiral)
library(dplyr)

input <- tribble(
  ~USUBJID, ~EXSTDTC, ~VISITNUM, ~EXTRT, ~EXLNKGRP, ~VISIT,
  "A001", "2015-01-10", 1, "DRUG A", "VAC 1", "VISIT 1",
  "A001", "2015-01-11", 2, "DRUG A", "VAC 2", "VISIT 2",
  "A001", "2015-01-12", 3, "DRUG B", "VAC 3", "VISIT 3",
  "A002", "2015-01-13", 1, "DRUG B", "VAC 1", "VISIT 1",
  "A002", "2015-01-14", 2, "DRUG C", "VAC 2", "VISIT 2"
)

adsl <- tribble(
  ~USUBJID, ~SEX, ~AGE,
  "A001", "MALE", 23,
  "A002", "FEMALE", 26,
)

derive_vars_vaxdt(
  dataset = input,
  dataset_adsl = adsl,
  by_vars = exprs(USUBJID, VISITNUM),
  order = exprs(USUBJID, VISITNUM, VISIT, EXSTDTC)
)
}
\seealso{
Other der_var: 
\code{\link{derive_var_aval_adis}()},
\code{\link{derive_vars_crit}()},
\code{\link{derive_vars_event_flag}()},
\code{\link{derive_vars_max_flag}()},
\code{\link{derive_vars_merged_vaccine}()},
\code{\link{derive_vars_params}()}
}
\author{
Vikram S
}
\concept{der_var}
\keyword{der_var}
