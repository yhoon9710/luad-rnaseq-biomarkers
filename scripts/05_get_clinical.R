# Download TCGA-LUAD patient-level survival data directly from the GDC API.
# Why: the clinical columns attached to the expression object are missing follow-up
# time for most living patients (which silently drops them from survival analysis),
# and TCGAbiolinks::GDCquery_clinic() currently errors on this project.
suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(tidyverse)
})

fields <- c("submitter_id",
            "demographic.vital_status", "demographic.days_to_death", "demographic.age_at_index",
            "diagnoses.days_to_last_follow_up", "diagnoses.days_to_death",
            "diagnoses.ajcc_pathologic_stage", "diagnoses.age_at_diagnosis",
            "follow_ups.days_to_follow_up")

resp <- POST("https://api.gdc.cancer.gov/cases",
             body = list(filters = list(op = "=",
                                        content = list(field = "project.project_id",
                                                       value = "TCGA-LUAD")),
                         fields = paste(fields, collapse = ","),
                         format = "JSON", size = "2000"),
             encode = "json")
stop_for_status(resp)
hits <- fromJSON(content(resp, "text", encoding = "UTF-8"),
                 simplifyVector = FALSE)$data$hits
message("Cases returned by GDC: ", length(hits))

num   <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))
maxna <- function(v) { v <- v[!is.na(v)]; if (length(v)) max(v) else NA_real_ }
first <- function(v) { v <- v[!is.na(v)]; if (length(v)) v[[1]] else NA }

clin <- map_dfr(hits, function(h) {
  dem   <- h$demographic %||% list()
  diags <- h$diagnoses   %||% list()
  fus   <- h$follow_ups  %||% list()
  tibble(
    patient          = h$submitter_id,
    vital_status     = dem$vital_status %||% NA_character_,
    days_to_death    = maxna(c(num(dem$days_to_death),
                               map_dbl(diags, ~ num(.x$days_to_death)))),
    days_last_contact = maxna(c(map_dbl(diags, ~ num(.x$days_to_last_follow_up)),
                                map_dbl(fus,   ~ num(.x$days_to_follow_up)))),
    age              = first(c(num(dem$age_at_index),
                               map_dbl(diags, ~ num(.x$age_at_diagnosis)) / 365.25)),
    stage_raw        = first(map_chr(diags, ~ .x$ajcc_pathologic_stage %||% NA_character_))
  )
}) %>%
  mutate(
    event = as.integer(vital_status == "Dead"),
    time_months = ifelse(event == 1, days_to_death, days_last_contact) / 30.44,
    # order matters: "Stage I" is a prefix of II/III/IV
    stage = case_when(str_detect(stage_raw, "Stage IV")  ~ "IV",
                      str_detect(stage_raw, "Stage III") ~ "III",
                      str_detect(stage_raw, "Stage II")  ~ "II",
                      str_detect(stage_raw, "Stage I")   ~ "I",
                      TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(event), !is.na(time_months), time_months > 0) %>%
  distinct(patient, .keep_all = TRUE)

dir.create("data", showWarnings = FALSE)
write_tsv(clin, "data/clinical_luad.tsv")

message("Patients with usable survival data: ", nrow(clin),
        " | alive: ", sum(clin$event == 0), " | dead: ", sum(clin$event == 1),
        " | median follow-up (months): ", round(median(clin$time_months), 1))
