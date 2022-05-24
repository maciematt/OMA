
library(dplyr)


## Pulling in the data for microarray meta-analysis testing.

d <- GEOquery::getGEO("GSE32924")[[1]]

d_expr <- d@assayData$exprs %>% as.data.frame %>% tibble::rownames_to_column(var = "ID") %>% OMA::translate_probes(probe_id = "ID")

d_info <- d@phenoData@data %>%
  transmute(
    sample_name = geo_accession,
    patient_id = characteristics_ch1.1 %>% stringr::str_remove(".*Patient "), sample_source = "skin", tissue = "skin",
    subject_id = characteristics_ch1.1 %>% stringr::str_remove(".*Patient "),
    sampling_time = "No Info",
    sample_id = geo_accession,
    sample_pathology = case_when(
      source_name_ch1 == "Normal" ~ "normal control",
      source_name_ch1 == "AL" ~ "lesional",
      source_name_ch1 == "ANL" ~ "non-lesional"
    ),
    disease_state = ifelse(stringr::str_detect(source_name_ch1, "Normal"), "normal control", "atopic dermatitis")
  )


## saveRDS(list(d_expr, d_info), file = "../inst/testdata/gse32924.rds")
saveRDS(list(d_expr, d_info), file = "testdata/gse32924.rds")

