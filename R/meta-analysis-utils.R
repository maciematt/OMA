
#' @import dplyr
#' @export translate_probes
#' @export nsfilter_oma



translate_probes <- function (d_expr, probe_id, from_id = "PROBEID", to_id = "ENTREZID", target_name = NULL, weighted_average = FALSE, annotation_db = hgu133plus2.db::hgu133plus2.db) {

  #' @param d_expr - A dataframe of expression where the first column contains the probe ids.
  #' @param probe_id - The name of the column that holds the probe ids.
  #' @param to_id - id type to which probes will be translated.
  #' @param target_name - the name of the column with the translated ID's (by default same as `to_id`)

  if (is.null(target_name))
    target_name <- to_id

  p2e <- AnnotationDbi::mapIds(annotation_db, keys = d_expr %>% pull(!!probe_id), column = to_id, keytype = from_id, multiVals = "list")
  p2e <- 1:length(p2e) %>% lapply(function (x) {
    tibble(probe = names(p2e)[x], targetid = p2e[[x]], inv_weight = length(p2e[[x]]))
  }) %>% bind_rows %>% tidyr::drop_na()

  if (isTRUE(weighted_average)) {
    d_expr_target <- p2e %>% left_join(d_expr %>% rename(probe = !!rlang::sym(probe_id))) %>% mutate(across(-c(probe, targetid, inv_weight), ~ .x * 1/inv_weight)) %>% select(-probe, -inv_weight) %>% group_by(targetid) %>% summarize(across(everything(), sum)) %>% left_join(
      p2e %>% tidyr::drop_na() %>% select(-probe) %>% group_by(targetid) %>% summarize(sum_weight = sum(1/inv_weight)), by = "targetid"
    ) %>% mutate(across(-c(targetid, sum_weight), ~ .x / sum_weight)) %>% rename(!!quo_name(target_name) := targetid) %>% select(-sum_weight)
  } else {
    d_expr_target <- p2e %>% left_join(d_expr %>% rename(probe = !!rlang::sym(probe_id))) %>% select(-probe, -inv_weight) %>% rename(!!quo_name(target_name) := targetid)
  }

  return(d_expr_target)

}


nsfilter_oma <- function (data_expr, data_pheno, samples, annotation) {

  d_bes <- Biobase::ExpressionSet(data_expr[, samples] %>% as.matrix, Biobase::AnnotatedDataFrame(as.data.frame(data_pheno[samples, ])), annotation = annotation)

  d_bes_f <- genefilter::nsFilter(
    d_bes, require.entrez = TRUE, remove.dupEntrez = TRUE,
    var.func = IQR, var.filter = TRUE, var.cutoff = 0.1, filterByQuantile = TRUE,
    feature.exclude = "^AFFX"
  )

  return (as.matrix(exprs(d_bes_f$eset)))
}
