## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
library(admiral)
library(admiralvaccine)
library(admiraldev)
library(pharmaversesdtm)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(stringr)
library(tidyr)
library(tibble)

data("face_vaccine")
data("suppface_vaccine")
data("ex_vaccine")
data("suppex_vaccine")
data("vs_vaccine")
data("admiralvaccine_adsl")

face <- convert_blanks_to_na(face_vaccine)
ex <- convert_blanks_to_na(ex_vaccine)
vs <- convert_blanks_to_na(vs_vaccine)
suppface <- convert_blanks_to_na(suppface_vaccine)
suppex <- convert_blanks_to_na(suppex_vaccine)
adsl <- convert_blanks_to_na(admiralvaccine_adsl)

## ----eval=TRUE----------------------------------------------------------------
face <- face %>%
  filter(FACAT == "REACTOGENICITY" & grepl("ADMIN|SYS", FASCAT)) %>%
  mutate(FAOBJ = str_to_upper(FAOBJ)) %>%
  metatools::combine_supp(suppface)
ex <- metatools::combine_supp(ex, suppex)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  face,
  display_vars = exprs(USUBJID, FAOBJ, FATESTCD, FACAT, FASCAT, FATPTREF, FADTC)
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_vars_merged_vaccine(
  dataset = face,
  dataset_ex = ex,
  by_vars_sys = exprs(USUBJID, FATPTREF = EXLNKGRP),
  by_vars_adms = exprs(USUBJID, FATPTREF = EXLNKGRP, FALOC = EXLOC, FALAT = EXLAT),
  ex_vars = exprs(EXTRT, EXDOSE, EXSEQ, EXSTDTC, EXENDTC, VISIT, VISITNUM)
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, FATESTCD, FATPTREF, EXTRT)
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RFSTDTC, RFENDTC)

adface <- derive_vars_merged(
  face,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, RFSTDTC, RFENDTC)
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_fever_records(
  dataset = adface,
  dataset_source = ungroup(vs),
  filter_source = VSCAT == "REACTOGENICITY" & VSTESTCD == "TEMP",
  faobj = "FEVER"
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, FATESTCD, FAORRES, VSSTRESN),
  filter = FAOBJ == "FEVER"
)

## ----eval=TRUE----------------------------------------------------------------
adface <- adface %>%
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = FADTC
  ) %>%
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = FADTC,
    highest_imputation = "n"
  )

## ----eval=TRUE----------------------------------------------------------------
adface <- adface %>%
  mutate(RFSTDTC = as.Date(RFSTDTC)) %>%
  derive_vars_dy(reference_date = RFSTDTC, source_vars = exprs(ADT))

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FATPTREF, FAOBJ, ADT, ADTM, ADY)
)

## ----eval=TRUE----------------------------------------------------------------
period_ref <- create_period_dataset(
  dataset = adsl,
  new_vars = exprs(
    APERSDT = APxxSDT, APEREDT = APxxEDT, TRTA = TRTxxA,
    TRTP = TRTxxP
  )
)

adface <- derive_vars_joined(
  adface,
  dataset_add = period_ref,
  by_vars = get_admiral_option("subject_keys"),
  filter_join = ADT >= APERSDT & ADT <= APEREDT,
  join_type = "all"
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, APERSDT, APEREDT, TRTA, TRTP)
)

## ----eval=TRUE----------------------------------------------------------------
sev_to_numeric <- function(x, y) {
  case_when(
    x == "NONE" ~ 0,
    x == "MILD" ~ 1,
    x == "MODERATE" ~ 2,
    x == "SEVERE" ~ 3,
    TRUE ~ y
  )
}

## ----eval=TRUE----------------------------------------------------------------
adface <- adface %>%
  mutate(
    AVALC = as.character(FASTRESC),
    AVAL = suppressWarnings(as.numeric(FASTRESN)),
    AVAL = sev_to_numeric(AVALC, AVAL),
    ATPTREF = FATPTREF,
    ATPT = FATPT,
    ATPTN = FATPTNUM
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, AVAL, AVALC, ATPTREF, ATPTN)
)

## ----eval = TRUE--------------------------------------------------------------
adface <- adface %>% derive_var_extreme_flag(
  by = exprs(STUDYID, USUBJID, FATPTREF, FAOBJ, FATESTCD, FATPTNUM),
  order = exprs(STUDYID, USUBJID, FATPTREF, FAOBJ, FATESTCD, FATPTNUM, FAEVAL),
  new_var = ANL01FL,
  mode = "first",
  true_value = "Y",
  false_value = NA_character_
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, ATPTREF, FATESTCD, FATEST, AVAL, ANL01FL),
  filter = ANL01FL == "Y"
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_diam_to_sev_records(
  dataset = adface,
  filter_add = ANL01FL == "Y",
  diam_code = "DIAMETER",
  faobj_values = c("REDNESS", "SWELLING"),
  testcd_sev = "SEV",
  test_sev = "Severity/Intensity",
  none = 0,
  mild = 2,
  mod = 5,
  sev = 10
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, ATPTREF, FATESTCD, FATEST, AVAL, AVALC),
  filter = FATESTCD == "SEV"
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_extreme_records(
  dataset = adface,
  dataset_add = adface,
  filter_add = FATESTCD == "SEV" & ANL01FL == "Y",
  by_vars = exprs(USUBJID, FAOBJ, ATPTREF),
  order = exprs(AVAL),
  check_type = "none",
  mode = "last",
  set_values_to = exprs(
    FATEST = "Maximum Severity",
    FATESTCD = "MAXSEV"
  )
)

adface <- derive_extreme_records(
  dataset = adface,
  dataset_add = adface,
  filter_add = FAOBJ %in% c("REDNESS", "SWELLING") & FATESTCD == "DIAMETER" & ANL01FL == "Y",
  by_vars = exprs(USUBJID, FAOBJ, FALNKGRP),
  order = exprs(AVAL),
  check_type = "none",
  mode = "last",
  set_values_to = exprs(
    FATEST = "Maximum Diameter",
    FATESTCD = "MAXDIAM"
  )
)

adface <- derive_extreme_records(
  dataset = adface,
  dataset_add = adface,
  filter_add = FAOBJ == "FEVER" & ANL01FL == "Y",
  by_vars = exprs(USUBJID, FAOBJ, ATPTREF),
  order = exprs(VSSTRESN),
  check_type = "none",
  mode = "last",
  set_values_to = exprs(
    FATEST = "Maximum Temperature",
    FATESTCD = "MAXTEMP"
  )
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, ATPTREF, FATESTCD, FATEST, AVAL, AVALC),
  filter = FATESTCD %in% c("MAXSEV", "MAXDIAM", "MAXTEMP")
)

## ----eval=TRUE, include=FALSE-------------------------------------------------
lookup_dataset <- tribble(
  ~FATESTCD, ~PARAMCD, ~PARAMN, ~FATEST, ~FAOBJ,
  "SEV", "SEVREDN", 1, "Severity/Intensity", "REDNESS",
  "DIAMETER", "DIARE", 2, "Diameter", "REDNESS",
  "MAXDIAM", "MDIRE", 3, "Maximum Diameter", "REDNESS",
  "MAXTEMP", "MAXTEMP", 4, "Maximum Temperature", "FEVER",
  "OCCUR", "OCFEVER", 5, "Occurrence Indicator", "FEVER",
  "OCCUR", "OCERYTH", 6, "Occurrence Indicator", "ERYTHEMA",
  "MAXSEV", "MAXSWEL", 7, "Maximum Severity", "SWELLING",
  "MAXSEV", "MAXREDN", 8, "Maximum Severity", "REDNESS",
  "MAXSEV", "MAXSFAT", 9, "Maximum Severity", "FATIGUE",
  "MAXSEV", "MAXSHEA", 10, "Maximum Severity", "HEADACHE",
  "MAXSEV", "MSEVNWJP", 11, "Maximum Severity", "NEW OR WORSENED JOINT PAIN",
  "MAXSEV", "MSEVNWMP", 12, "Maximum Severity", "NEW OR WORSENED MUSCLE PAIN",
  "OCCUR", "OCISR", 13, "Occurrence Indicator", "REDNESS",
  "OCCUR", "OCINS", 14, "Occurrence Indicator", "SWELLING",
  "OCCUR", "OCPIS", 15, "Occurrence Indicator", "PAIN AT INJECTION SITE",
  "OCCUR", "OCFATIG", 16, "Occurrence Indicator", "FATIGUE",
  "OCCUR", "OCHEAD", 17, "Occurrence Indicator", "HEADACHE",
  "OCCUR", "OCCHILLS", 18, "Occurrence Indicator", "CHILLS",
  "OCCUR", "OCDIAR", 19, "Occurrence Indicator", "DIARRHEA",
  "OCCUR", "OCCNWJP", 20, "Occurrence Indicator", "NEW OR WORSENED JOINT PAIN",
  "OCCUR", "OCCNWMP", 21, "Occurrence Indicator", "NEW OR WORSENED MUSCLE PAIN",
  "SEV", "SEVSWEL", 22, "Severity/Intensity", "SWELLING",
  "SEV", "SEVPIS", 23, "Severity/Intensity", "PAIN AT INJECTION SITE",
  "SEV", "SEVFAT", 24, "Severity/Intensity", "FATIGUE",
  "SEV", "SEVHEAD", 25, "Severity/Intensity", "HEADACHE",
  "SEV", "SEVDIAR", 26, "Severity/Intensity", "DIARRHEA",
  "SEV", "SEVNWJP", 27, "Severity/Intensity", "NEW OR WORSENED JOINT PAIN",
  "SEV", "SEVNWMP", 28, "Severity/Intensity", "NEW OR WORSENED MUSCLE PAIN",
  "MAXDIAM", "MDISW", 29, "Maximum Diameter", "SWELLING",
  "MAXSEV", "MAXSPIS", 30, "Maximum Severity", "PAIN AT INJECTION SITE",
  "OCCUR", "OCCVOM", 31, "Occurrence Indicator", "VOMITING",
  "DIAMETER", "DIASWEL", 32, "Diameter", "SWELLING"
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_vars_params(
  dataset = adface,
  lookup_dataset = lookup_dataset,
  merge_vars = exprs(PARAMCD, PARAMN)
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, ATPTREF, FATEST, PARAMCD, PARAM, PARAMN, PARCAT1, PARCAT2)
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_vars_max_flag(
  dataset = adface,
  flag1 = "ANL02FL",
  flag2 = "ANL03FL"
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, ATPTREF, AVAL, AVALC, ANL01FL, ANL02FL)
)

## ----eval=TRUE----------------------------------------------------------------
adface <- derive_vars_event_flag(
  dataset = adface,
  by_vars = exprs(USUBJID, FAOBJ, ATPTREF),
  aval_cutoff = 2.5,
  new_var1 = EVENTFL,
  new_var2 = EVENTDFL
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, FAOBJ, ATPTREF, AVAL, AVALC, EVENTFL, EVENTDFL)
)

## ----eval=TRUE----------------------------------------------------------------
adface <- post_process_reacto(
  dataset = adface,
  filter_dataset = FATESTCD %in% c("MAXDIAM", "MAXSEV", "MAXTEMP") |
    (FATESTCD %in% c("OCCUR", "SEV") & FAOBJ %in% c("FEVER", "REDNESS", "SWELLING"))
)

## ----eval=TRUE----------------------------------------------------------------
adsl <- adsl %>%
  convert_blanks_to_na() %>%
  filter(!is.na(USUBJID))

adface <- derive_vars_merged(
  dataset = adface,
  dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
  by_vars = get_admiral_option("subject_keys")
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adface,
  display_vars = exprs(USUBJID, TRTSDT, TRTEDT, AGE, SEX)
)

