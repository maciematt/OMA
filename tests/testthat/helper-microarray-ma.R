
library(dplyr)


## Here reading in a pre-stored GEO dataset, preparing it for analysis, then
## proceeding to run the analysis, first by hand, then by OMA. This should pro-
## duce the same (exremely similar) results for the test to work correctly.


d <- system.file("testdata", "gse32924.rds", package = "OMA") %>% readRDS

d_expr <- d[[1]]
d_info <- d[[2]]
rownames(d_info) <- d_info$sample_name
d <- NULL

keep_samples <- intersect(d_info %>% filter(sample_pathology != "normal control") %>% pull(sample_name), d_expr %>% colnames)

d_expr_t <- d_expr %>% ungroup %>% select(1, d_info %>% pull(sample_id)) %>% select(-ENTREZID) %>% as.matrix %>% t
colnames(d_expr_t) <- d_expr$ENTREZID


d_expr_ready <- d_expr %>% select(-1) %>% as.data.frame
rownames(d_expr_ready) <- d_expr$ENTREZID

## At this point the data is ready; now limma by hand

mm <- model.matrix(~ 0 + d_info[keep_samples, "sample_pathology", drop = TRUE])
colnames(mm) <- levels(d_info[keep_samples, "sample_pathology", drop = TRUE] %>% factor) %>% make.names


con <- limma::makeContrasts(contrasts = "lesional - non.lesional", levels = mm)

corfit <- limma::duplicateCorrelation(d_expr_ready[, keep_samples], mm, block = d_info[keep_samples, "patient_id", drop = TRUE])
lmfit_realized <- limma::lmFit(d_expr_ready[, keep_samples], design = mm, block = d_info[keep_samples, "patient_id", drop = TRUE], correlation = corfit$result$consensus.correlation)


limma_manual <- limma::eBayes(
  limma::contrasts.fit(
    lmfit_realized,
    con
  )
)

## Done with limma by hand; time for OMA

d_ready <- prepare_for_dge(list(list(sample_data = d_info[keep_samples, ], gene_expression_data = d_expr_t[keep_samples, ])), "GSE10761", contrast = list(variable = "sample_pathology", active = "lesional", reference = "non-lesional"), blocks = "subject_id", scrutinize_blocks = TRUE)

d_dge <- run_dge(d_ready)

## Both done, now can run tests
