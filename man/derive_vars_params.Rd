% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_params.R
\name{derive_vars_params}
\alias{derive_vars_params}
\title{Assigning Parameter Variables}
\usage{
derive_vars_params(dataset, lookup_dataset, merge_vars)
}
\arguments{
\item{dataset}{Input dataset
Input dataset is expected to have variables \code{USUBJID},\code{FAOBJ},
\code{FACAT}, \code{FATESTCD} and \code{FATEST}}

\item{lookup_dataset}{lookup dataset containing \code{PARAMCD} values for every
unique \code{FATESTCD} and \code{FAOBJ}
lookup dataset is expected to have variables \code{FATEST}, \code{PARAMCD},
\code{FATESTCD}, \code{FAOBJ} and one entry for every unique
\code{FATESTCD} and \code{FAOBJ}}

\item{merge_vars}{List of Variables need to be merged from lookup dataset}
}
\value{
The output dataset contains all observations and variables of the
input dataset along with \code{PARAM},\code{PARAMCD},\code{PARCAT1},\code{PARCAT2},\code{PARAMN}
}
\description{
Creating \code{PARAMCD} from lookup dataset and assigning \code{PARAM},\code{PARAMN},\code{PARCAT1},
\code{PARCAT2} variables
}
\details{
A lookup dataset is required with \code{PARAMCD} values for every combination
of \code{FATEST} & \code{FAOBJ}.
\code{PARAMCD} \code{PARAMN} \code{PARAMN} \code{PARCAT1} \code{PARCAT2} values can be assigned
from lookup dataset.

\if{html}{\out{<div class="sourceCode">}}\preformatted{ if `PARAMN` not assigned in lookup dataset then
 `PARAMN` is assigned with a unique number for every unique PARAM value.
 if `PARAM` value not assigned in lookup dataset then
 `PARAM` value is a combination of `FAOBJ` `FATEST` `FASTRESU` `FALOC`
 `FADIR` `FALAT`
 if `PARCAT1` value not assigned in lookup dataset then
 `PARCAT1` is assigned as `FACAT`
 if `PARCAT2` value not assigned in lookup dataset then
 `PARCAT2` is assigned as `FASCAT`
}\if{html}{\out{</div>}}
}
\examples{

library(admiral)
library(tibble)
library(dplyr)

lookup_dataset <- tibble::tribble(
  ~FATESTCD, ~PARAMCD, ~PARAMN, ~FATEST, ~FAOBJ,
  "SEV", "SEVREDN", 1, "Severity", "Redness",
  "DIAMETER", "DIARE", 2, "Diameter", "Redness",
  "MAXDIAM", "MDIRE", 3, "Maximum Diameter cm", "Redness",
  "MAXTEMP", "MAXTEMP", 4, "Maximum Temperature", "Fever",
  "OCCUR", "OCFEVER", 5, "Occurrence Indicator", "Fever",
  "OCCUR", "OCERYTH", 6, "Occurrence Indicator", "Erythema",
  "SEV", "SEVPAIN", 7, "Severity", "Pain at Injection site",
  "OCCUR", "OCPAIN", 8, "Occurrence Indicator", "Pain at Injection site",
  "OCCUR", "OCSWEL", 9, "Occurrence Indicator", "Swelling"
)

input <- tibble::tribble(
  ~USUBJID, ~FACAT, ~FASCAT, ~FATESTCD, ~FAOBJ, ~FATEST, ~FALOC, ~FALAT,
  "ABC101", "REACTO", "ADMIN", "SEV", "Redness", "Severity", "ARM", "LEFT",
  "ABC101", "REACTO", "ADMIN", "DIAMETER", "Redness", "Diameter", "ARM", "RIGHT",
  "ABC101", "REACTO", "ADMIN", "MAXDIAM", "Redness", "Maximum Diameter", NA, NA,
  "ABC101", "REACTO", "SYSTEMIC", "MAXTEMP", "Fever", "Maximum Temp", NA, NA,
  "ABC101", "REACTO", "SYSTEMIC", "OCCUR", "Fever", "Occurrence", NA, NA,
  "ABC101", "REACTO", "ADMIN", "OCCUR", "Erythema", "Occurrence", NA, NA,
  "ABC101", "REACTO", "ADMIN", "SEV", "Swelling", "Severity", NA, NA,
  "ABC101", "REACTO", "ADMIN", "OCCUR", "Swelling", "Occurrence", NA, NA,
  "ABC101", "REACTO", "ADMIN", "OCCUR", "Swelling", "Occurrence", NA, NA
)

derive_vars_params(
  dataset = input,
  lookup_dataset = lookup_dataset,
  merge_vars = exprs(PARAMCD, PARAMN)
)

}
\seealso{
Other der_var: 
\code{\link{derive_var_aval_adis}()},
\code{\link{derive_vars_crit}()},
\code{\link{derive_vars_event_flag}()},
\code{\link{derive_vars_max_flag}()},
\code{\link{derive_vars_merged_vaccine}()},
\code{\link{derive_vars_vaxdt}()}
}
\author{
Dhivya Kanagaraj
}
\concept{der_var}
\keyword{der_var}
