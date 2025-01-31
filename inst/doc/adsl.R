## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(admiral)
library(admiralvaccine)
library(pharmaversesdtm)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(stringr)
library(admiraldev)

data("dm_vaccine")
data("ex_vaccine")

dm <- convert_blanks_to_na(dm_vaccine)
ex <- convert_blanks_to_na(ex_vaccine)

## ----eval=TRUE----------------------------------------------------------------
adsl <- dm %>%
  select(-DOMAIN)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, RFSTDTC, COUNTRY, AGE, SEX, RACE, ETHNIC, ARM, ACTARM)
)

## ----eval=TRUE----------------------------------------------------------------
adsl <- dm %>%
  mutate(
    TRT01P = substring(ARM, 1, 9),
    TRT02P = substring(ARM, 11, 100)
  ) %>%
  derive_vars_merged(
    dataset_add = ex,
    filter_add = EXLNKGRP == "VACCINATION 1",
    new_vars = exprs(TRT01A = EXTRT),
    by_vars = get_admiral_option("subject_keys")
  ) %>%
  derive_vars_merged(
    dataset_add = ex,
    filter_add = EXLNKGRP == "VACCINATION 2",
    new_vars = exprs(TRT02A = EXTRT),
    by_vars = get_admiral_option("subject_keys")
  )

## ----eval=TRUE----------------------------------------------------------------
# impute start and end time of exposure to first and last respectively, do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN"
  )
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "VACCINE"))) &
      !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = get_admiral_option("subject_keys")
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "VACCINE"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = get_admiral_option("subject_keys")
  )

## ----eval=TRUE----------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))

## ----eval=TRUE----------------------------------------------------------------
adsl <- adsl %>%
  derive_var_trtdurd()

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, RFSTDTC, TRTSDTM, TRTSDT, TRTEDTM, TRTEDT, TRTDURD)
)

## ----eval=TRUE----------------------------------------------------------------
adsl <- derive_var_merged_exist_flag(
  dataset = adsl,
  dataset_add = ex,
  by_vars = exprs(STUDYID, USUBJID),
  new_var = SAFFL,
  condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "VACCINE")))
) %>%
  mutate(
    PPROTFL = "Y"
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTSDT, ARM, ACTARM, SAFFL, PPROTFL)
)

## ----eval=TRUE----------------------------------------------------------------
adsl <- derive_vars_vaxdt(
  dataset = ex,
  dataset_adsl = adsl,
  by_vars = exprs(USUBJID, VISITNUM),
  order = exprs(USUBJID, VISITNUM, VISIT, EXSTDTC)
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, VAX01DT, VAX02DT)
)

## ----eval=TRUE----------------------------------------------------------------
adsl <- adsl %>%
  mutate(
    AP01SDT = VAX01DT,
    AP01EDT = if_else(!is.na(VAX02DT), VAX02DT - 1, as.Date(RFPENDTC)),
    AP02SDT = if_else(!is.na(VAX02DT), VAX02DT, NA_Date_),
    AP02EDT = if_else(!is.na(AP02SDT), as.Date(RFPENDTC), NA_Date_)
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, AP01SDT, AP01EDT, AP02SDT, AP02EDT)
)

