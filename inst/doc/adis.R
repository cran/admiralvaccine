## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(admiral)
link <- function(text, url) {
  return(
    paste0(
      "[", text, "]",
      "(", url, ")"
    )
  )
}
dyn_link <- function(text,
                     base_url,
                     relative_url = "",
                     # Change to TRUE when admiral adopts multiversion docs
                     is_multiversion = FALSE,
                     multiversion_default_ref = "main") {
  url <- paste(base_url, relative_url, sep = "/")
  if (is_multiversion) {
    url <- paste(
      base_url,
      Sys.getenv("BRANCH_NAME", multiversion_default_ref),
      relative_url,
      sep = "/"
    )
  }
  return(link(text, url))
}
# Other variables
admiral_homepage <- "https://pharmaverse.github.io/admiral"

## ----eval=TRUE, message=FALSE, warning=FALSE----------------------------------
library(admiral)
library(dplyr)
library(lubridate)
library(admiraldev)
library(admiralvaccine)
library(pharmaversesdtm)
library(metatools)

# Load source datasets
data("is_vaccine")
data("suppis_vaccine")
data("admiralvaccine_adsl")

# Convert blanks into NA
is <- convert_blanks_to_na(is_vaccine)
suppis <- convert_blanks_to_na(suppis_vaccine)
adsl <- convert_blanks_to_na(admiralvaccine_adsl)

## ----eval=TRUE----------------------------------------------------------------
is_suppis <- metatools::combine_supp(is, suppis)

## ----eval=TRUE----------------------------------------------------------------
adis <- is_suppis %>%
  mutate(
    AVISITN = as.numeric(VISITNUM),
    AVISIT = case_when(
      VISITNUM == 10 ~ "Visit 1",
      VISITNUM == 20 ~ "Visit 2",
      VISITNUM == 30 ~ "Visit 3",
      VISITNUM == 40 ~ "Visit 4",
      is.na(VISITNUM) ~ NA_character_
    ),
    ATPTN = as.numeric(VISITNUM / 10),
    ATPT = case_when(
      VISITNUM == 10 ~ "Visit 1 (Day 1)",
      VISITNUM == 20 ~ "Visit 2 (Day 31)",
      VISITNUM == 30 ~ "Visit 3 (Day 61)",
      VISITNUM == 40 ~ "Visit 4 (Day 121)",
      is.na(VISITNUM) ~ NA_character_
    ),
    ATPTREF = case_when(
      VISITNUM %in% c(10, 20) ~ "FIRST TREATMENT",
      VISITNUM %in% c(30, 40) ~ "SECOND TREATMENT",
      is.na(VISITNUM) ~ NA_character_
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, AVISIT, AVISITN, ATPT, ATPTN, ATPTREF)
)

## ----eval=TRUE----------------------------------------------------------------
# ADT derivation
# Add also PPROTFL from ADSL (to avoid additional merges) in order to derive
# PPSRFL at step 11.
adis <- derive_vars_dt(
  dataset = adis,
  new_vars_prefix = "A",
  dtc = ISDTC,
  highest_imputation = "M",
  date_imputation = "mid",
  flag_imputation = "none"
) %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(RFSTDTC, PPROTFL),
    by_vars = get_admiral_option("subject_keys")
  ) %>%
  mutate(
    ADT = as.Date(ADT),
    RFSTDTC = as.Date(RFSTDTC)
  ) %>%
  # ADY derivation
  derive_vars_dy(
    reference_date = RFSTDTC,
    source_vars = exprs(ADT)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, ISDTC, RFSTDTC, ADT, ADY, PPROTFL)
)

## ----eval=TRUE----------------------------------------------------------------
# Create record duplication in order to plot both original and LOG10 parameter values.
# Add also records related to 4fold.
# Please, keep or modify PARAM values according to your purposes.

is_log <- adis %>%
  mutate(
    DERIVED = "LOG10",
    ISSEQ = NA_real_
  )

is_4fold <- adis %>%
  mutate(
    DERIVED = "4FOLD",
    ISSEQ = NA_real_
  )

is_log_4fold <- adis %>%
  mutate(
    DERIVED = "LOG10 4FOLD",
    ISSEQ = NA_real_
  )

adis <- bind_rows(adis, is_log, is_4fold, is_log_4fold) %>%
  arrange(STUDYID, USUBJID, !is.na(DERIVED), ISSEQ) %>%
  mutate(DERIVED = if_else(is.na(DERIVED), "ORIG", DERIVED))


adis <- adis %>%
  mutate(
    # PARAMCD: for log values, concatenation of L and ISTESTCD.
    PARAMCD = case_when(
      DERIVED == "ORIG" ~ ISTESTCD,
      DERIVED == "LOG10" ~ paste0(ISTESTCD, "L"),
      DERIVED == "4FOLD" ~ paste0(ISTESTCD, "F"),
      # As per CDISC rule, PARAMCD should be 8 characters long. Please, adapt if needed
      DERIVED == "LOG10 4FOLD" ~ paste0(substr(ISTESTCD, 1, 6), "LF")
    )
  )


# Update param_lookup dataset with your PARAM values.
param_lookup <- tribble(
  ~PARAMCD, ~PARAM, ~PARAMN,
  "J0033VN", "J0033VN Antibody", 1,
  "I0019NT", "I0019NT Antibody", 2,
  "M0019LN", "M0019LN Antibody", 3,
  "R0003MA", "R0003MA Antibody", 4,
  "J0033VNL", "LOG10 (J0033VN Antibody)", 11,
  "I0019NTL", "LOG10 (I0019NT Antibody)", 12,
  "M0019LNL", "LOG10 (M0019LN Antibody)", 13,
  "R0003MAL", "LOG10 (R0003MA Antibody)", 14,
  "J0033VNF", "4FOLD (J0033VN Antibody)", 21,
  "I0019NTF", "4FOLD (I0019NT Antibody)", 22,
  "M0019LNF", "4FOLD (M0019LN Antibody)", 23,
  "R0003MAF", "4FOLD (R0003MA Antibody)", 24,
  "J0033VLF", "LOG10 4FOLD (J0033VN Antibody)", 31,
  "I0019NLF", "LOG10 4FOLD (I0019NT Antibody)", 32,
  "M0019LLF", "LOG10 4FOLD (M0019LN Antibody)", 33,
  "R0003MLF", "LOG10 4FOLD (R0003MA Antibody)", 34
)

adis <- derive_vars_merged_lookup(
  dataset = adis,
  dataset_add = param_lookup,
  new_vars = exprs(PARAM, PARAMN),
  by_vars = exprs(PARAMCD)
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, PARAMCD, PARAM, PARAMN)
)

## ----eval=TRUE----------------------------------------------------------------
adis <- adis %>%
  mutate(
    PARCAT1 = ISCAT,
    # Please, define your additional cutoff values. Delete if not needed.
    CUTOFF02 = 4,
    CUTOFF03 = 8
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, PARCAT1, CUTOFF02, CUTOFF03)
)

## ----eval=TRUE----------------------------------------------------------------
adis_or <- adis %>%
  filter(DERIVED == "ORIG") %>%
  derive_var_aval_adis(
    lower_rule = ISLLOQ / 2,
    middle_rule = ISSTRESN,
    upper_rule = ISULOQ,
    round = 2
  )

adis_log_or <- adis %>%
  filter(DERIVED == "LOG10") %>%
  derive_var_aval_adis(
    lower_rule = log10(ISLLOQ / 2),
    middle_rule = log10(ISSTRESN),
    upper_rule = log10(ISULOQ),
    round = 2
  )

adis_4fold <- adis %>%
  filter(DERIVED == "4FOLD") %>%
  derive_var_aval_adis(
    lower_rule = ISLLOQ,
    middle_rule = ISSTRESN,
    upper_rule = ISULOQ,
    round = 2
  )

adis_log_4fold <- adis %>%
  filter(DERIVED == "LOG10 4FOLD") %>%
  derive_var_aval_adis(
    lower_rule = log10(ISLLOQ),
    middle_rule = log10(ISSTRESN),
    upper_rule = log10(ISULOQ),
    round = 2
  )

adis <- bind_rows(adis_or, adis_log_or, adis_4fold, adis_log_4fold) %>%
  mutate(
    # AVALU derivation (please delete if not needed for your study)
    AVALU = ISSTRESU,

    # SERCAT1 derivation
    SERCAT1 = case_when(
      ISBLFL == "Y" & !is.na(AVAL) & !is.na(ISLLOQ) & AVAL < ISLLOQ ~ "S-",
      ISBLFL == "Y" & !is.na(AVAL) & !is.na(ISLLOQ) & AVAL >= ISLLOQ ~ "S+",
      ISBLFL == "Y" & (is.na(AVAL) | is.na(ISLLOQ)) ~ "UNKNOWN"
    )
  )


# Update param_lookup2 dataset with your SERCAT1N values.
param_lookup2 <- tribble(
  ~SERCAT1, ~SERCAT1N,
  "S-", 1,
  "S+", 2,
  "UNKNOWN", 3,
  NA_character_, NA_real_
)

adis <- derive_vars_merged_lookup(
  dataset = adis,
  dataset_add = param_lookup2,
  new_vars = exprs(SERCAT1N),
  by_vars = exprs(SERCAT1)
)


# DTYPE derivation.
# Please update code when <,<=,>,>= are present in your lab results (in ISSTRESC)

if (any(names(adis) == "ISULOQ") == TRUE) {
  adis <- adis %>%
    mutate(DTYPE = case_when(
      DERIVED %in% c("ORIG", "LOG10") & !is.na(ISLLOQ) &
        ((ISSTRESN < ISLLOQ) | grepl("<", ISORRES)) ~ "HALFLLOQ",
      DERIVED %in% c("ORIG", "LOG10") & !is.na(ISULOQ) &
        ((ISSTRESN > ISULOQ) | grepl(">", ISORRES)) ~ "ULOQ",
      TRUE ~ NA_character_
    ))
}

if (any(names(adis) == "ISULOQ") == FALSE) {
  adis <- adis %>%
    mutate(DTYPE = case_when(
      DERIVED %in% c("ORIG", "LOG10") & !is.na(ISLLOQ) &
        ((ISSTRESN < ISLLOQ) | grepl("<", ISORRES)) ~ "HALFLLOQ",
      TRUE ~ NA_character_
    ))
}

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, AVAL, AVALU, DTYPE, SERCAT1, SERCAT1N)
)

## ----eval=TRUE----------------------------------------------------------------
# ABLFL derivation
adis <- restrict_derivation(
  adis,
  derivation = derive_var_extreme_flag,
  args = params(
    by_vars = exprs(STUDYID, USUBJID, PARAMN),
    order = exprs(STUDYID, USUBJID, VISITNUM, PARAMN),
    new_var = ABLFL,
    mode = "first"
  ),
  filter = VISITNUM == 10
) %>%
  # BASE derivation
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMN),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  # BASETYPE derivation
  derive_basetype_records(
    basetypes = exprs("VISIT 1" = AVISITN %in% c(10, 30))
  ) %>%
  arrange(STUDYID, USUBJID, !is.na(DERIVED), ISSEQ)


# BASECAT derivation
adis <- adis %>%
  mutate(
    BASECAT1 = case_when(
      !grepl("L", PARAMCD) & BASE < 10 ~ "Titer value < 1:10",
      !grepl("L", PARAMCD) & BASE >= 10 ~ "Titer value >= 1:10",
      grepl("L", PARAMCD) & BASE < 10 ~ "Titer value < 1:10",
      grepl("L", PARAMCD) & BASE >= 10 ~ "Titer value >= 1:10"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, ABLFL, BASE, BASETYPE, BASECAT1)
)

## ----eval=TRUE----------------------------------------------------------------
adis <- restrict_derivation(adis,
  derivation = derive_var_chg,
  filter = AVISITN > 10
) %>%
  restrict_derivation(
    derivation = derive_var_analysis_ratio,
    args = params(
      numer_var = AVAL,
      denom_var = BASE
    ),
    filter = AVISITN > 10
  ) %>%
  arrange(STUDYID, USUBJID, DERIVED, ISSEQ)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, CHG, R2BASE)
)

## ----eval=TRUE----------------------------------------------------------------
# STEP 9 Derivation of CRITyFL and CRITyFN ----
adis <- restrict_derivation(
  dataset = adis,
  derivation = derive_vars_crit_flag,
  args = params(
    crit_nr = 1,
    condition = AVAL >= ISLLOQ,
    description = "Titer >= ISLLOQ",
    values_yn = TRUE,
    create_numeric_flag = TRUE
  ),
  filter = !is.na(AVAL)
) %>%
  arrange(STUDYID, USUBJID, DERIVED, ISSEQ)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, CRIT1, CRIT1FL, CRIT1FN)
)

## ----eval=TRUE----------------------------------------------------------------
period_ref <- create_period_dataset(
  dataset = adsl,
  new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT, TRTA = TRTxxA, TRTP = TRTxxP)
)

adis <- derive_vars_joined(
  adis,
  dataset_add = period_ref,
  by_vars = get_admiral_option("subject_keys"),
  filter_join = ADT >= APERSDT & ADT <= APEREDT,
  join_type = "all"
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, TRTP, TRTA)
)

## ----eval=TRUE----------------------------------------------------------------
adis <- adis %>%
  mutate(PPSRFL = if_else(VISITNUM %in% c(10, 30) & PPROTFL == "Y", "Y", NA_character_))

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, TRTP, TRTA)
)

## ----eval=TRUE----------------------------------------------------------------
# Get list of ADSL variables not to be added to ADIS
vx_adsl_vars <- exprs(RFSTDTC, PPROTFL)

adis <- derive_vars_merged(
  dataset = adis,
  dataset_add = select(adsl, !!!negate_vars(vx_adsl_vars)),
  by_vars = get_admiral_option("subject_keys")
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adis,
  display_vars = exprs(USUBJID, VISITNUM, ISTEST, ISORRES, AGE, COUNTRY, ARM, ACTARM)
)

