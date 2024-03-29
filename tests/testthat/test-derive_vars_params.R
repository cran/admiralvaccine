## Test 1: Check if the variables in the lookup dataset getting merged properly

test_that("derive_vars_params Test 1: Check if the variables in the lookup dataset getting merged
           properly", {
  lookup_dataset <- tibble::tribble(
    ~FATESTCD,    ~PARAMCD,   ~FAOBJ,
    "SEV",        "SEVREDN",  "Redness",
    "DIAMETER",   "DIARE",    "Redness",
    "MAXDIAM",    "MDIRE",    "Redness",
    "MAXTEMP",    "MAXTEMP",  "Fever",
    "OCCUR",      "OCFEVER",  "Fever",
    "OCCUR",      "OCERYTH",  "Erythema",
    "SEV",        "SEVPAIN",  "Pain at Injection site",
    "OCCUR",      "OCPAIN",   "Pain at Injection site",
    "OCCUR",      "OCSWEL",   "Swelling"
  )

  input <- tibble::tribble(
    ~USUBJID, ~FACAT, ~FASCAT, ~FATESTCD, ~FAOBJ, ~FATEST,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "SEV", "Redness", "Severity",
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "DIAMETER", "Redness", "Diameter",
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "MAXDIAM", "Redness", "Maximum Diameter",
    "ABC101", "REACTOGENICITY", "SYSTEMIC", "MAXTEMP", "Fever", "Maximum Temp",
    "ABC101", "REACTOGENICITY", "SYSTEMIC", "OCCUR", "Fever", "Occurrence",
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "OCCUR", "Erythema", "Occurrence",
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "SEV", "Swelling", "Severity",
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "OCCUR", "Swelling", "Occurrence"
  )

  expout1 <- input %>%
    left_join(lookup_dataset,
      by = c("FATESTCD", "FAOBJ")
    ) %>%
    convert_na_to_blanks() %>%
    mutate(
      PARCAT1 = FACAT,
      PARCAT2 = FASCAT,
      PARAM = gsub(
        "(^[[:space:]]+|[[:space:]]+$)", "",
        str_to_sentence(paste0(FAOBJ, " ", FATEST))
      )
    )
  expout2 <- expout1 %>%
    distinct(PARAM, .keep_all = FALSE) %>%
    mutate(PARAMN = seq_len(n()))

  expected_output <- merge(expout1, expout2, by = "PARAM", all.x = TRUE) %>%
    convert_blanks_to_na()

  actual_output <- derive_vars_params(
    dataset = input,
    lookup_dataset = lookup_dataset,
    merge_vars = exprs(PARAMCD)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("USUBJID", "PARAM", "PARAMCD", "PARCAT1", "PARCAT2", "PARAMN")
  )
})

## Test 2: Checking whether PARAM  getting concatenated with only the existed variables

test_that("derive_vars_params Test 2: Checking whether PARAM  getting concatenated with
                    only the existed variables", {
  lookup_dataset <- tibble::tribble(
    ~FATESTCD,    ~PARAMCD,   ~FAOBJ,
    "SEV",        "SEVREDN",  "Redness",
    "DIAMETER",   "DIARE",    "Redness",
    "MAXDIAM",    "MDIRE",    "Redness",
    "MAXTEMP",    "MAXTEMP",  "Fever",
    "OCCUR",      "OCFEVER",  "Fever",
    "OCCUR",      "OCERYTH",  "Erythema",
    "SEV",        "SEVPAIN",  "Pain at Injection site",
    "OCCUR",      "OCPAIN",   "Pain at Injection site",
    "OCCUR",      "OCSWEL",   "Swelling"
  )

  input <- tibble::tribble(
    ~USUBJID, ~FACAT, ~FASCAT, ~FATESTCD, ~FAOBJ, ~FATEST, ~FALAT, ~FALOC, ~FADIR,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "SEV", "Redness", "Severity",
    "DELTOID MUSCLE", "LEFT", NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "DIAMETER", "Redness", "Diameter",
    "DELTOID MUSCLE", "LEFT", NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "MAXDIAM", "Redness", "Maximum Diameter",
    "DELTOID MUSCLE", "RIGHT", NA,
    "ABC101", "REACTOGENICITY", "SYSTEMIC", "MAXTEMP", "Fever", "Maximum Temp",
    NA, NA, NA,
    "ABC101", "REACTOGENICITY", "SYSTEMIC", "OCCUR", "Fever", "Occurrence",
    NA, NA, NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "OCCUR", "Erythema", "Occurrence",
    "RIGHT", NA, NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "SEV", "Swelling", "Severity", NA,
    NA, NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "OCCUR", "Swelling", "Occurrence",
    NA, "RIGHT", NA
  )



  expout1 <- input %>%
    left_join(lookup_dataset,
      by = c("FATESTCD", "FAOBJ")
    ) %>%
    mutate(
      PARCAT1 = FACAT,
      PARCAT2 = FASCAT,
      PARAM = ""
    ) %>%
    unite(PARAM, FAOBJ, FATEST, FADIR, FALOC, FALAT,
      sep = " ",
      na.rm = TRUE, remove = FALSE
    ) %>%
    mutate(PARAM = str_to_sentence(PARAM))



  expout2 <- expout1 %>%
    distinct(PARAM, .keep_all = FALSE) %>%
    mutate(PARAMN = seq_len(n()))

  expected_output <- merge(expout1, expout2, by = "PARAM", all.x = TRUE)

  actual_output <- derive_vars_params(
    dataset = as.data.frame(input),
    lookup_dataset = lookup_dataset,
    merge_vars = exprs(PARAMCD)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("USUBJID", "PARAM", "PARAMCD", "PARCAT1", "PARCAT2", "PARAMN")
  )
})


## Test 3: Checking whether PARAM, PARCAT1 and PARCAT2 getting concatenated with
## only the existed variables

test_that("derive_vars_params Test 3: Checking whether PARAM, PARCAT1 and PARCAT2 getting
           concatenated with only the existed variables", {
  lookup_dataset <- tibble::tribble(
    ~FATESTCD, ~PARAMCD, ~FAOBJ, ~PARAMN,
    "SEV", "SEVREDN", "Redness", 1,
    "DIAMETER", "DIARE", "Redness", 2,
    "MAXDIAM", "MDIRE", "Redness", 3,
    "MAXTEMP", "MAXTEMP", "Fever", 4,
    "OCCUR", "OCFEVER", "Fever", 5,
    "OCCUR", "OCERYTH", "Erythema", 6,
    "SEV", "SEVPAIN", "Pain at Injection site", 7,
    "OCCUR", "OCPAIN", "Pain at Injection site", 8,
    "OCCUR", "OCSWEL", "Swelling", 9
  )

  input <- tibble::tribble(
    ~USUBJID, ~FACAT, ~FASCAT, ~FATESTCD, ~FAOBJ, ~FATEST, ~FALAT, ~FALOC, ~FADIR,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "SEV", "Redness", "Severity",
    "DELTOID MUSCLE", "LEFT", NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "DIAMETER", "Redness", "Diameter",
    "DELTOID MUSCLE", "LEFT", NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "MAXDIAM", "Redness", "Maximum Diameter",
    "DELTOID MUSCLE", "RIGHT", NA,
    "ABC101", "REACTOGENICITY", "SYSTEMIC", "MAXTEMP", "Fever", "Maximum Temp",
    NA, NA, NA,
    "ABC101", "REACTOGENICITY", "SYSTEMIC", "OCCUR", "Fever", "Occurrence",
    NA, NA, NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "OCCUR", "Erythema", "Occurrence",
    "RIGHT", NA, NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "SEV", "Swelling", "Severity", NA,
    NA, NA,
    "ABC101", "REACTOGENICITY", "ADMIN-SITE", "OCCUR", "Swelling", "Occurrence",
    NA, "RIGHT", NA
  )

  expected_output <- input %>%
    left_join(lookup_dataset,
      by = c("FATESTCD", "FAOBJ")
    ) %>%
    mutate(
      PARCAT1 = FACAT,
      PARCAT2 = FASCAT,
      PARAM = ""
    ) %>%
    unite(PARAM, FAOBJ, FATEST, FADIR, FALOC, FALAT,
      sep = " ",
      na.rm = TRUE, remove = FALSE
    ) %>%
    mutate(PARAM = str_to_sentence(PARAM))

  actual_output <- derive_vars_params(
    dataset = as.data.frame(input),
    lookup_dataset = lookup_dataset,
    merge_vars = exprs(PARAMCD, PARAMN)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("USUBJID", "PARAM", "PARAMCD", "PARCAT1", "PARCAT2", "PARAMN")
  )
})
