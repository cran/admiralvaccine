% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_var_aval_adis.R
\name{derive_var_aval_adis}
\alias{derive_var_aval_adis}
\title{Derive AVAL variable for ADIS ADaM domain}
\usage{
derive_var_aval_adis(
  dataset,
  lower_rule,
  middle_rule,
  upper_rule = NULL,
  round = NULL
)
}
\arguments{
\item{dataset}{Input dataset.}

\item{lower_rule}{Derivation rule when \code{ISSTRESN} value is below \code{ISLLOQ}.
When \code{ISSTRESN} is missing, the inequality in \code{ISORRES} is checked for the derivation.}

\item{middle_rule}{Derivation rule when \code{ISSTRESN} value is greater than \code{ISLLOQ} and
lower than \code{ISULOQ}.
If \code{ISULOQ} is not present, derivation rule when \code{ISSTRESN} is greater than \code{ISLLOQ}.
When \code{ISSTRESN} is missing, the inequality in \code{ISORRES} is checked for the derivation.}

\item{upper_rule}{Derivation rule when \code{ISSTRESN} value is greater than \code{ISULOQ}.
This is an optional argument since \code{ISULOQ} may not be present.
When \code{ISSTRESN} is missing, the inequality in \code{ISORRES}. is checked for the derivation.
Default value is \code{NULL}.}

\item{round}{Rounding for \code{AVAL} variable. An integer argument which specifies
the number of decimals displayed.
Default value is \code{NULL}.}
}
\value{
Dataset with \code{AVAL} variable derived.
}
\description{
Derive \code{AVAL} variable for Laboratory Immunology Data ADaM domain.
A common rule has been decided for its derivation, based on \code{ISLLOQ}, \code{ISULOQ} and \code{ISORRES}
when both \code{ISLLOQ} and \code{ISULOQ} are present.
If \code{ISULOQ} is not present, the variables used are \code{ISLLOQ} and \code{ISORRES}.
Please, refers to arguments description for additional details.
}
\examples{
library(tibble)
library(admiral)
library(admiraldev)
library(dplyr)
library(rlang)

input <- tribble(
  ~USUBJID, ~AVISITN, ~PARAMCD, ~PARAM, ~ISORRES, ~ISSTRESN, ~ISLLOQ, ~ISULOQ,
  "ABC-1001", 10, "J0033VN", "J0033VN Antibody", NA, NA, 2, 100,
  "ABC-1001", 10, "I0019NT", "I0019NT Antibody", "3", 3.0, 4, 200,
  "ABC-1001", 10, "M0019LN", "M0019LN Antibody", ">150", NA, 8, 150,
  "ABC-1001", 10, "R0003MA", "R0003MA Antibody", "140.5", 140.5, 4, 120,
  "ABC-1001", 30, "J0033VN", "J0033VN Antibody", "2", 2.0, 2, 100,
  "ABC-1001", 30, "I0019NT", "I0019NT Antibody", NA, NA, 4, 200,
  "ABC-1001", 30, "M0019LN", "M0019LN Antibody", NA, NA, 8, 150,
  "ABC-1001", 30, "R0003MA", "R0003MA Antibody", "98.2", 98.2, 4, 120,
  "ABC-1001", 10, "J0033VNL", "LOG10 (J0033VN Antibody)", NA, NA, 2, 100,
  "ABC-1001", 10, "I0019NTL", "LOG10 (I0019NT Antibody)", "3", 3.0, 4, 200,
  "ABC-1001", 10, "M0019LNL", "LOG10 (M0019LN Antibody)", ">150", NA, 8, 150,
  "ABC-1001", 10, "R0003MAL", "LOG10 (R0003MA Antibody)", "140.5", 140.5, 4, 120,
  "ABC-1001", 30, "J0033VNL", "LOG10 (J0033VN Antibody)", "2", 2.0, 2, 100,
  "ABC-1001", 30, "I0019NTL", "LOG10 (I0019NT Antibody)", NA, NA, 4, 200,
  "ABC-1001", 30, "M0019LNL", "LOG10 (M0019LN Antibody)", NA, NA, 8, 150,
  "ABC-1001", 30, "R0003MAL", "LOG10 (R0003MA Antibody)", "98.2", 98.2, 4, 120,
  "ABC-1002", 10, "J0033VN", "J0033VN Antibody", "3", 3.0, 2, 100,
  "ABC-1002", 10, "I0019NT", "I0019NT Antibody", NA, NA, 4, 200,
  "ABC-1002", 10, "M0019LN", "M0019LN Antibody", NA, NA, 8, 150,
  "ABC-1002", 10, "R0003MA", "R0003MA Antibody", "48.9", 48.9, 4, 120,
  "ABC-1002", 30, "J0033VN", "J0033VN Antibody", NA, NA, 2, 100,
  "ABC-1002", 30, "I0019NT", "I0019NT Antibody", NA, NA, 4, 200,
  "ABC-1002", 30, "M0019LN", "M0019LN Antibody", "5", 5.0, 8, 150,
  "ABC-1002", 30, "R0003MA", "R0003MA Antibody", "228.1", 228.1, 4, 120
)

derive_var_aval_adis(
  dataset = input,
  lower_rule = ISLLOQ / 2,
  middle_rule = ISSTRESN,
  upper_rule = ISULOQ,
  round = 2
)

}
\seealso{
Other der_var: 
\code{\link{derive_vars_crit}()},
\code{\link{derive_vars_event_flag}()},
\code{\link{derive_vars_max_flag}()},
\code{\link{derive_vars_merged_vaccine}()},
\code{\link{derive_vars_params}()},
\code{\link{derive_vars_vaxdt}()}
}
\author{
Federico Baratin
}
\concept{der_var}
\keyword{der_var}
