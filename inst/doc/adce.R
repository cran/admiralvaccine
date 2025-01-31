## ----setup, include = FALSE---------------------------------------------------
# Other variables
admiral_homepage <- "https://pharmaverse.github.io/admiral"

## ----eval=TRUE, message=FALSE, warning=FALSE----------------------------------
library(admiraldev)
library(admiral)
library(dplyr)
library(lubridate)
library(admiralvaccine)
library(pharmaversesdtm)

data("ce_vaccine")
data("admiralvaccine_adsl")

adsl <- admiralvaccine_adsl
ce <- ce_vaccine

ce <- convert_blanks_to_na(ce)
adsl <- convert_blanks_to_na(adsl)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTSDT, TRTEDT, TRT01A, AP01SDT, AP01EDT, AP02SDT, AP02EDT)
)

## ----eval=TRUE----------------------------------------------------------------
adce <- ce %>%
  filter(CECAT == "REACTOGENICITY")

## ----eval=TRUE----------------------------------------------------------------
adsl2 <- adsl %>%
  select(-c(starts_with("AP") & ends_with("DTM")))

adperiods <- create_period_dataset(
  adsl2,
  new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT)
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adperiods,
  display_vars = exprs(USUBJID, APERIOD, APERSDT, APEREDT)
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(TRTSDT, TRTEDT)

adce <- adce %>%
  # join ADSL to CE
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = get_admiral_option("subject_keys")
  ) %>%
  derive_vars_dt(
    dtc = CESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "n"
  ) %>%
  derive_vars_dt(
    dtc = CEENDTC,
    new_vars_prefix = "AEN",
    highest_imputation = "n"
  ) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ASTDT, AENDT)
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adce,
  display_vars = exprs(USUBJID, TRTSDT, CESTDTC, CEENDTC, ASTDT, AENDT, ASTDY, AENDY)
)

## ----eval=TRUE----------------------------------------------------------------
adce <-
  derive_vars_joined(
    adce,
    dataset_add = adperiods,
    by_vars = get_admiral_option("subject_keys"),
    filter_join = ASTDT >= APERSDT & ASTDT <= APEREDT,
    join_type = "all"
  ) %>%
  mutate(
    APERSTDY = as.integer(ASTDT - APERSDT) + 1,
    AREL = CEREL
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adce,
  display_vars = exprs(USUBJID, TRTSDT, ASTDT, AENDT, ASTDY, AENDY, APERIOD, APERSDT, APERSTDY)
)

## ----eval=TRUE----------------------------------------------------------------
adce <- adce %>%
  mutate(
    ASEV = CESEV,
    ASEVN = as.integer(factor(ASEV,
      levels = c("MILD", "MODERATE", "SEVERE", "DEATH THREATENING")
    ))
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, APERIOD),
      order = exprs(desc(ASEVN), ASTDY, CEDECOD),
      new_var = AOCC01FL,
      mode = "first"
    ),
    filter = !is.na(APERIOD) & !is.na(ASEV)
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adce,
  display_vars = exprs(
    USUBJID, TRTSDT, ASTDT, APERIOD, APERSDT,
    APERSTDY, CEDECOD, ASEVN, AOCC01FL, CESEQ
  )
)

## ----eval=TRUE----------------------------------------------------------------
adce <- adce %>%
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(CEDECOD, CELAT, CETPTREF, APERIOD),
    check_type = "error"
  ) %>%
  derive_vars_duration(
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT,
    in_unit = "days",
    out_unit = "days",
    add_one = TRUE,
    trunc_out = FALSE
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adce,
  display_vars = exprs(
    USUBJID, TRTSDT, ASTDT, APERIOD, APERSDT,
    APERSTDY, CEDECOD, ASEVN, AOCC01FL, CESEQ, ASEQ
  )
)

## ----eval=TRUE----------------------------------------------------------------
adsl_list <- adsl %>%
  select(STUDYID, USUBJID, TRT01A, TRT01P, AGE, AGEU, SEX, RACE, COUNTRY, ETHNIC, SITEID, SUBJID)

adce <- adce %>%
  derive_vars_merged(
    dataset_add = adsl_list,
    by_vars = get_admiral_option("subject_keys")
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adce,
  display_vars = exprs(
    USUBJID, TRTSDT, ASTDT, APERIOD, APERSDT, APERSTDY,
    CEDECOD, ASEVN, AOCC01FL, CESEQ, ASEQ, AGE, SEX
  )
)

