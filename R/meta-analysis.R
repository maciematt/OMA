

#' @import dplyr
#' @export prepare_for_dge
#' @export run_dge
#' @export run_dgsva
#' @export run_ma
#' @export



prepare_for_dge <- function (
  dataset,
  dataset_name,
  contrast,
  contrast_type = c("binary", "numeric"),
  technology = c("microarray", "RNA-seq"),
  sample_col = "sample_name",
  covariates = NULL,
  scrutinize_covariates = TRUE,
  blocks = NULL,
  scrutinize_blocks = FALSE,
  sample_filters = NULL,
  use_sva = FALSE,
  sva_method = c("leek"),
  min_genes_to_keep_dataset = 5,
  min_samples_per_level = 2,
  min_samples = 5,
  fix_NA = c("none", "mean", "median"),
  filter_genes = TRUE,
  annotation_db = hgu133plus2.db::hgu133plus2.db,
  translate_to_entrez = TRUE
) {

  #' Remember to pre-merge any synonymous contrast levels before using OMA!
  #'
  #' dataset - this needs to ba a list representing a dataset object with two
  #'   fields: "gene_expression_data", and "sample_data".
  #' dataset_name - this needs to ba a vector of dataset names - same order as
  #'   datasets above.
  #' contrast - if `contrast_type` is "binary" a three-member list, if 
  #'   "numeric" one-memebered list, where (in both cases) the first field is 
  #'   the name of the  variable to use from "pheno", then in "binary" the
  #'   second is the name of the active level in the contrast, and the third is
  #'   the name of the reference level.
  #' sample_filters - a quosure of filters; previously this was a string that
  #'   would be passed into dplyr's `filter_()`; however, `filter_()` was de-
  #'   precated, and this was quite error-prone anyway, as it required a double
  #'   qotations in and similar risky business in the case of strings. So now a
  #'   full expression is expected that goes directly into dplyr's `filter()` -
  #'   send it in inside `rlang::quo()` which will make it a quosure.
  #' response_variable - informs which variable is the one to use as a contrast


  technology <- match.arg(technology)
  contrast_type <- match.arg(contrast_type)
  sva_method <- match.arg(sva_method)


  if (!is.null(sample_filters) & !rlang::is_quosure(sample_filters)) { stop("`sample_filters` needs to be a quosure - create it using the `create_filter` function!") }

  if (contrast_type == "binary")
    if (!is.list(contrast) | (length(contrast) != 3) | !all(names(contrast) == c("variable", "active", "reference"))) { stop("For a binary contrast, `contrast` needs to be a 3-value list, with the following properties: \"variable\" (the name of the contrast variable in the data), \"active\" (the active level(s) of the contrast), \"reference\" (reference level(s) of the contrast).") }
  else if (contrast_type == "numeric")
    if (!is.list(contrast) | (length(contrast) != 1) | !all(names(contrast) == c("variable"))) { stop("For a numeric contrast, `contrast` needs to be a 1-value list, with the following property: \"variable\" (the name of the contrast variable in the data).") }

  if (!is.null(covariates) & !is.vector(covariates)) { stop("`covariates` needs to be a vector of covariates that will be used in dataset!") }
  if (!is.null(blocks) & !is.vector(blocks)) { stop("`blocks` needs to be a vector of blocks that will be used across all datasets!") }
  if(is.vector(blocks) & (length(blocks) != 1)) { stop("If `blocks` is a vector, than it needs to be a vector of length 1 (only a single variable can be used to identify duplicates).") }

  if (!is.null(annotation_db) & attr(class(annotation_db), "package") != "AnnotationDbi") { stop("`annotation_db` needs to be an AnnotationDbi-packaged object!") }

  if (!all(names(dataset) %in% c("gene_expression_data", "sample_data"))) { stop(paste0("each dataset needs to be a list with two dataframes, one for gene expression (a named list element called \"gene_expression_data\"), and sample data (another named list element called \"sample_data\"). this is not the case in ", dataset_name)) }
  if (!is.data.frame(dataset$gene_expression_data) & !is.matrix(dataset$gene_expression_data)) { stop(paste0("a \"gene_expression_data\" field is present in ", dataset_name, ", but it's not a data frame or matrix.")) }
  if (!is.data.frame(dataset$sample_data)) { stop(paste0("a \"sample_data\" field is present in ", dataset_name, ", but it's not a data frame.")) }
  if (!(sample_col %in% colnames(dataset$sample_data))) { stop(paste0("the \"sample_data\" dataframe needs to have a \"", sample_col, "\" column corresponding to the rownames in the expression data. this column wasn't found in ", dataset_name, "!")) }


  ## testing
  # contrast <- list(
  #   variable = "disease_state",
  #   active = "atopic dermatitis",
  #   reference = "normal control"
  # )

  response_variable <- contrast[["variable"]]

  if (contrast_type == "binary") {
    active_levels <- contrast[["active"]]
    reference_levels <- contrast[["reference"]]
  }


  fix_NA <- match.arg(fix_NA)


  ## testing
  # covariates <- "gender"

  if (!is.null(covariates) & scrutinize_covariates) {
    remaining_covariates <- intersect(covariates, colnames(dataset$sample_data))
    if (length(remaining_covariates) > 0)
      if (contrast_type == "binary")
        remaining_covariates <- remaining_covariates[remaining_covariates %>% sapply(function (x) dataset$sample_data %>% filter(!!rlang::sym(response_variable) %in% c(active_levels, reference_levels)) %>% pull(x) %>% as.character %>% na.omit %>% unique %>% length >= 2)]
      else
        ## for a binary contrast still requesting that a given covariate is represented in at least two observations
        remaining_covariates <- remaining_covariates[remaining_covariates %>% sapply(function (x) dataset$sample_data %>% tidyr::drop_na(!!rlang::sym(response_variable)) %>% pull(x) %>% as.character %>% na.omit %>% unique %>% length >= 2)]

    if (!identical(sort(covariates), sort(remaining_covariates)))
      stop(paste0("Some of the requested covariates cannot be included.\nRequested covariates:\n", paste0(covariates, collapse = ", "), "\nCovariates with >2 levels:\n", paste0(remaining_covariates, collapse = ", ")))
  }


  if (!is.null(blocks) & scrutinize_blocks) {
    remaining_blocks <- intersect(blocks, colnames(dataset$sample_data))
    if (length(remaining_blocks) > 0)
      if (contrast_type == "binary")
        remaining_blocks <- remaining_blocks[remaining_blocks %>% sapply(function (x) dataset$sample_data %>% filter(!!rlang::sym(response_variable) %in% c(active_levels, reference_levels)) %>% pull(x) %>% as.character %>% unique %>% length >= 2)]
      else
        remaining_blocks <- remaining_blocks[remaining_blocks %>% sapply(function (x) dataset$sample_data %>% filter(dplyr::drop_na(!!rlang::sym(response_variable))) %>% pull(x) %>% length >= 2)]
    if (length(remaining_blocks) == 0)
      remaining_blocks <- NULL
    blocks <- remaining_blocks
  }


  info <- list(
    run_date = date(),
    min_genes_to_keep_dataset = min_genes_to_keep_dataset,
    min_samples_per_level = min_samples_per_level,
    min_samples = min_samples,
    contrast = contrast,
    contrast_type = contrast_type,
    technology = technology,
    scrutinize_covariates = scrutinize_covariates,
    covariates = covariates,
    scrutinize_blocks = scrutinize_blocks,
    blocks = blocks,
    fix_NA = fix_NA,
    dataset_name = dataset_name,
    filter = filter,
    annotation_db = annotation_db,
    use_sva = use_sva,
    sva_method = sva_method,
    translate_to_entrez = translate_to_entrez
  )


  if (!all(covariates %in% colnames(dataset$sample_data))) { stop(paste0("Covariates specified for dataset ", dataset_name, " that don't exist in its pheno data")) }
  if (!all(blocks %in% colnames(dataset$sample_data))) { stop(paste0("`blocks` specified for dataset ", dataset_name, " that don't exist in its pheno data")) }

  cat(paste0("dataset: ", dataset_name, "\n"))


  data_obj <- list()
  data_obj$status <- list(
    NA_present = FALSE,
    expr_dims_eq_0 = c(FALSE, FALSE),
    keep = TRUE
  )
  data_obj$pheno <- dataset$sample_data
  data_obj$covariates <- covariates
  data_obj$block <- blocks

  print("data_obj$status:")
  print(data_obj$status)


  ## Checking if pheno meets the criteria ----------------------------------- ##

  if (!is.null(sample_filters))
    data_obj[["pheno"]] <- data_obj[["pheno"]] %>% filter(!!sample_filters) %>% as.data.frame

  if (contrast_type == "binary")
    data_obj[["pheno"]] <- data_obj[["pheno"]] %>% mutate(!!response_variable := factor(!!rlang::sym(response_variable)))
  else if (contrast_type == "numeric")
    data_obj[["pheno"]] <- data_obj[["pheno"]] %>% mutate(!!response_variable := as.numeric(!!rlang::sym(response_variable)))

  ## For uniformity, changing the key sample column name to "sample_name"
  data_obj[["pheno"]] <- data_obj[["pheno"]] %>% rename(sample_name := !!sample_col)

  c_n <- unname(unlist(data_obj[["pheno"]]$sample_name))
  rownames(data_obj[["pheno"]]) <- c_n

  PHENO_EXPR_SAMPLE_NAME_INTERSECT <- intersect(c_n, dataset$gene_expression_data %>% rownames)

  # data_obj[["pheno"]] <- data_obj[["pheno"]][PHENO_EXPR_SAMPLE_NAME_INTERSECT, ]
  data_obj[["pheno"]] <- data_obj[["pheno"]] %>%
    filter(sample_name %in% PHENO_EXPR_SAMPLE_NAME_INTERSECT) %>%
    arrange(match(sample_name, PHENO_EXPR_SAMPLE_NAME_INTERSECT)) %>%
    as.data.frame
  
  rownames(data_obj[["pheno"]]) <- PHENO_EXPR_SAMPLE_NAME_INTERSECT
  
  print(levels(as.factor(unname(unlist(data_obj[["pheno"]][, response_variable])))))

  print("response variable bit")
  # print(data_obj[["pheno"]] %>% filter(!!rlang::sym(response_variable) %in% c(active_levels, reference_levels)) %>% pull(!!response_variable) %>% as.character %>% unique)

  data_obj$status$less_than_2_response_levels <- ifelse(
    contrast_type == "numeric",
    FALSE, # should be NA, going with FALSE to slot into the existing logic
    data_obj[["pheno"]] %>% filter(!!rlang::sym(response_variable) %in% c(active_levels, reference_levels)) %>% pull(!!response_variable) %>% as.character %>% unique %>% length < 2
  )
  print("less than 2 required levels?")
  print(data_obj$status$less_than_2_response_levels)

  data_obj$status$less_than_min_samples_per_level <- ifelse(
    contrast_type == "numeric",
    FALSE, # should be NA, going with FALSE to slot into the existing logic
    (function () {
      tmp_dat <- data_obj[["pheno"]] %>% filter(!!rlang::sym(response_variable) %in% c(active_levels, reference_levels)) %>% pull(!!response_variable) %>% as.character
      lvls <- tmp_dat %>% unique
      sapply(lvls, function (lvl) {
        length(tmp_dat[tmp_dat == lvl]) < min_samples_per_level
      }) %>% unlist %>% unname %>% any
    })()
  )

  data_obj$status$less_than_min_samples <- data_obj[["pheno"]] %>% filter(!is.na(!!rlang::sym(response_variable))) %>% pull(!!response_variable) %>% length < min_samples

  if (data_obj$status$less_than_2_response_levels) {
    cat(paste0("Less than 2 levels present in the response variable! Eliminating dataset ", dataset_name, "\n"))
  }
  if (data_obj$status$less_than_min_samples) {
    cat(paste0("Less than the minimal specified number (", min_samples, ") samples! Eliminating dataset ", dataset_name, "\n"))
  }
  if (data_obj$status$less_than_min_samples_per_level) {
    cat(paste0("Less than ", min_samples_per_level, " samples with a given response level present in the response variable! Eliminating dataset ", dataset_name, "\n"))
  }

  if (data_obj$status$less_than_2_response_levels | data_obj$status$less_than_min_samples_per_level | data_obj$status$less_than_min_samples) {
    data_obj$status$keep <- FALSE
    stop(paste0("Encountered a problem with the dataset:\n", data_obj$status))
  }


  c_n <- rownames(data_obj[["pheno"]])

  ## ------------------------------------------------------------------------ ##


  ## Expression data -------------------------------------------------------- ##

  data_obj$expr <- (dataset$gene_expression_data %>% as.data.frame)[PHENO_EXPR_SAMPLE_NAME_INTERSECT, ]


  if (any(is.na(data_obj[["expr"]])) | any(is.null(data_obj[["expr"]]))) {

    cat("\nNA values present...\n\n")
    data_obj$status$NA_present <- TRUE

    if (fix_NA == "mean") {

      gene_means <- lapply(data_obj[["expr"]], function (col) {
        if (any(is.na(col)) | any(is.null(col))) {
          col %>% mean(na.rm = TRUE)
	  ## col[!is.nan(col)] %>% unlist %>% unname %>% mean(na.rm = TRUE)
        } else { NULL }
      })
      names(gene_means) <- colnames(data_obj[["expr"]])
      gene_means[sapply(gene_means, is.null)] <- NULL
      gene_means <- unlist(gene_means)

      cat(paste0("NAs need fixing in ", length(gene_means), " genes\n"))
      # print(gene_means[1:10])

      cat("fixing NA's: ")
      for (gene_num in 1:length(gene_means)) {
        gene <- names(gene_means)[gene_num]
        if (gene_num %% 100 == 0)
          cat(paste0(gene_num, " "))
        data_obj[["expr"]][is.na(data_obj[["expr"]][, gene]), gene] <- gene_means[gene]
      }
      print("")

      cat(paste0("Checking if there are any NAs left... ", any(is.na(data_obj[["expr"]])), "\n"))
    }
  }


  r_n <- colnames(data_obj[["expr"]])


  ## c_n below comes from the names of samples in sample_data
  # data_obj[["expr"]] <- data_obj[["expr"]][c_n, ] %>% t %>% as.data.frame
  data_obj[["expr"]] <- data_obj[["expr"]] %>% t %>% as.data.frame

  # colnames(data_obj[["expr"]]) <- c_n
  colnames(data_obj[["expr"]]) <- PHENO_EXPR_SAMPLE_NAME_INTERSECT
  rownames(data_obj[["expr"]]) <- r_n

  ## EXPRESSION DATA FILTERING ---------------------------------------------- ##

  ## filter_genes ----------------------------------------------------------- ##
  if (isTRUE(filter_genes) & technology == "microarray") {
    ## nsFilter pollutes the session by attaching a lot of packages, so we call this garbage as a separate session
    annotation_txt <- annotation_db$packageName
    oma_tempdata <- file.path(tempdir(), "OMA-tempdata.rds")
    saveRDS(list(expr_data = data_obj[["expr"]], pheno_data = data_obj[["pheno"]], samples = PHENO_EXPR_SAMPLE_NAME_INTERSECT, annotation = annotation_txt), oma_tempdata)
    oma_tempscript <- file.path(tempdir(), "OMA-tempscript.R")
    nsfilter_run <- paste0('library(OMA)
    d <- readRDS("', oma_tempdata, '")
    out <- nsfilter_oma(d$expr_data, d$pheno_data, d$samples, d$annotation)
    saveRDS(out, "', oma_tempdata, '")
    ')
    writeLines(nsfilter_run, oma_tempscript)
    system(
      paste0("Rscript --vanilla --max-ppsize=500000 ", oma_tempscript)
    )
    data_obj[["expr"]] <- readRDS(oma_tempdata)
  }
  ## ------------------------------------------------------------------------ ##


  ## translate_to_entrez ---------------------------------------------------- ##
  if (isTRUE(translate_to_entrez)) {

    if (technology == "microarray") {
      d_ex <- bind_cols(tibble(probe = rownames(data_obj[["expr"]])), data_obj[["expr"]] %>% as.data.frame %>% as_tibble) %>% translate_probes(probe_id = "probe", annotation_db = annotation_db)
    } else if (technology == "RNA-seq") {
      d_ex <- bind_cols(tibble(gene = rownames(data_obj[["expr"]])), data_obj[["expr"]] %>% as.data.frame %>% as_tibble) %>% translate_probes(probe_id = "gene", from_id = "ENSEMBL", weighted_average = TRUE, annotation_db = org.Hs.eg.db::org.Hs.eg.db)
    }

    entrez_names <- d_ex %>% pull(1)
    data_obj[["expr"]] <- d_ex %>% select(-1) %>% as.data.frame %>% as.matrix
    rownames(data_obj[["expr"]]) <- entrez_names
  }
  ## ------------------------------------------------------------------------ ##

  ## ------------------------------------------------------------------------ ##



  # data_obj[["keys"]] <- rownames(data_obj[["expr"]])
  # names(data_obj[["keys"]]) <- rownames(data_obj[["expr"]])
  data_obj[["dataset_name"]] <- dataset_name


  c_n <- colnames(data_obj[["expr"]])
  r_n <- rownames(data_obj[["expr"]])
  temp <- data_obj[["expr"]] %>% t %>%
    as.data.frame %>% lapply(function (col) {
      as.numeric(col)
    }) %>% do.call(cbind, .) %>% as.matrix %>% t
  data_obj[["expr"]] <- temp
  colnames(data_obj[["expr"]]) <- c_n
  rownames(data_obj[["expr"]]) <- r_n


  ## ------------------------------------------------------------------------ ##


  print(paste0("dim: ", dim(data_obj[["expr"]])))

  if (0 %in% dim(data_obj[["expr"]])) {
    data_obj$status$expr_dims_eq_0 <- dim(data_obj[["expr"]]) == 0
    data_obj$status$keep <- FALSE
  }

  #   return(data_obj)
  # })


  if (isFALSE(data_obj$status$keep)) {
    warning(paste0("`status$keep` is FALSE for dataset ", dataset_name, "!\nHere's the full `status`:\n"))
    print(data_obj$status)
    return(NA)
  }


  ret_obj <- structure(
    list(
      info = info,
      data = data_obj,
      dataset_name = dataset_name,
      dataset_prep_info = data_obj$status
    ), class = "dge_ready_data"
  )

  return(ret_obj)

}



check_ge_obj <- function (ge_obj) {}



run_dge <- function (ge_object) {

  #' Here using limma for #1.
  #'
  #' @param `technology` needs to be a vector of "microarray" or "RNA-seq" values corresponding to the technologies in the corresponding members of `ge_object`. If left at NULL, all datasets in `ge_object` are assumed to be from microarrays.
  #'
  #' *Could be passing in a concrete contrast argument instead with a more
  #' complex contrast once needed (when/if implementing, remember that these
  #' will need to be supplied with `make.names`'d variable names) - for now
  #' working with simple contrasts, so don't need to implement this.
  #' *This expression should involve all levels from the "active_levels" and
  #' "reference_levels" lists.
  #'
  #' Remember to pre-merge any synonymous levels before using OMA!


  if (!is(ge_object, "dge_ready_data")) { stop("`ge_object` needs to be of type `dge_ready_data`!") }
  # lapply(ge_object, check_ge_obj)


  ge_info <- ge_object$info
  contrast_info <- ge_info$contrast

  if (ge_info$contrast_type == "binary") { ## this is the case only for binary contrasts
    active_nicename <- make.names(contrast_info[["active"]])
    reference_nicename <- make.names(contrast_info[["reference"]])
    contrast_nicename <- paste0(active_nicename, "-", reference_nicename)
  }


  # ge_results <- lapply(1:length(ge_object$data), function (d_num) {

  tech <- ge_info$technology

  ge_data <- ge_object$data


  ## Right up front, `model.matrix` is not able to automatically handle any rows with NA's.
  nonna_rows <- apply(ge_data$pheno[, c(contrast_info$variable, ge_data$covariates), drop = FALSE], 1, function(x) all(!is.na(x)))
  ge_data$expr <- ge_data$expr[, nonna_rows]
  ge_data$pheno <- ge_data$pheno[nonna_rows, ]


  ge_data$expr_ready <- ge_data$expr



  if (tech == "RNA-seq") {
    dlist <- edgeR::DGEList(counts = ge_data$expr_ready, group = ge_data$pheno[, contrast_info$variable, drop = TRUE])
    keep <- edgeR::filterByExpr(dlist, min.count = 0, min.total.count = 10*ncol(ge_data$expr_ready))
    dlist <- dlist[keep, , keep.lib.sizes = T]
    tmm <- edgeR::calcNormFactors(dlist, method = "TMM")
    ge_data$expr_ready <- tmm
  }

  ge_f <- ifelse(
    ge_info$contrast_type == "binary",
    paste0('~ 0 + ge_data[["pheno"]][, "', contrast_info$variable, '", drop = TRUE]'),
    paste0('~ ge_data[["pheno"]][, "', contrast_info$variable, '", drop = TRUE]')
  )
  if (!is.null(ge_data$covariates))
    ge_f <- paste0(ge_f, " + ", paste0(sapply(ge_data$covariates, function (x) paste0('ge_data$pheno[, "', x, '", drop = TRUE]')), collapse = " + "))
  ge_f <- as.formula(ge_f)

  # print(ge_data[["pheno"]][, contrast_info$variable])
  # print(contrast_info)
  # print(ge_data[["pheno"]][, "gender"])

  mm <- model.matrix(ge_f)

  if (ge_info$contrast_type == "binary") ## this is the case only for binary contrasts
    colnames(mm)[1:length(levels(ge_data[["pheno"]][, contrast_info$variable, drop = TRUE]))] <- levels(ge_data[["pheno"]][, contrast_info$variable, drop = TRUE])
  else if (ge_info$contrast_type == "numeric") ## for the binary contrast, we'll go with the simple variable name
    colnames(mm)[2] <- contrast_info$variable ## index=2 because of the interecept in the first position

  colnames(mm) <- make.names(colnames(mm))

  if (isTRUE(ge_info$use_sva)) {
    nonna_expr <- fix_missing_tx(ge_data$expr_ready, rownames(ge_data$expr_ready), colnames(ge_data$expr_ready))
    # n_sv <- sva::num.sv(ge_data$expr_ready, mm, method = ge_info$sva_method)
    n_sv <- sva::num.sv(nonna_expr, mm, method = ge_info$sva_method)
    if (n_sv == 0)
      stop("Can't run with automated SVA - number of surrogate variables picked by SVA is 0!")
    # sva_obj <- sva::sva(ge_data$expr_ready, mm, n.sv = n_sv)
    sva_obj <- sva::sva(nonna_expr, mm, n.sv = n_sv)
    sva_sv <- sva_obj$sv
    colnames(sva_sv) <- paste0("sv", seq(1, n_sv))
    mm <- cbind(mm, sva_sv)
  }

  print(colnames(mm))

  if (ge_info$contrast_type == "binary") ## this is the case only for binary contrasts
    con <- limma::makeContrasts(contrasts = contrast_nicename, levels = mm)


  if (tech == "RNA-seq")
    ge_data$expr_ready <- limma::voom(ge_data$expr_ready, mm)


  lmfit_partial <- purrr::safely(limma::lmFit) %>% purrr::partial(ge_data$expr_ready, design = mm)
  lmfit_realized <- lmfit_partial()$result

  if (!is.null(ge_data$block)) {
    corfit <- limma::duplicateCorrelation(ge_data$expr_ready, mm, block = ge_data$pheno[, ge_data$block, drop = TRUE])
    lmfit_partial_cor <- lmfit_partial %>% purrr::partial(block = ge_data$pheno[, ge_data$block, drop = TRUE], correlation = corfit$result$consensus.correlation)
    lmfit_partial_cor_realized <- lmfit_partial_cor()
    if (is.null(lmfit_partial_cor_realized$error)) {
      lmfit_realized <- lmfit_partial_cor_realized$result
    }
  }

  if (ge_info$contrast_type == "binary") ## this is the case only for binary contrasts
    limma_realized <- limma::contrasts.fit(
      lmfit_realized,
      con
    )
  else if (ge_info$contrast_type == "numeric")
    limma_realized <- lmfit_realized

  limma_out <- limma::eBayes(limma_realized)

  # limma::topTable(limma_out, number = Inf) %>% tibble::rownames_to_column(var = "gene") %>% head

  ## Note that this is focused on the contrast variable (all the `[, 1]`s or `[, 2]`s below)
  if (ge_info$contrast_type == "binary")
    ix_contrast <- 1
  else if (ge_info$contrast_type == "numeric")
    ix_contrast <- 2

  diff_data_out <- tibble(
    gene = rownames(limma_out$coefficients),
    coeff = limma_out$coefficients[, ix_contrast] %>% as.vector,
    p_value = limma_out$p.value[, ix_contrast] %>% as.vector,
    cohen_d = limma_out$coefficients[, ix_contrast] / sqrt(limma_out$s2.post), # As recommended by Gordon Smyth: https://support.bioconductor.org/p/71747/
    variance = limma_out$s2.post * limma_out$stdev.unscaled[, ix_contrast]^2 # --||-- in https://support.bioconductor.org/p/70175/
  )

  if (ge_info$contrast_type == "binary") { ## this is the case only for binary contrasts
    diff_data_out$active_n <- ge_data$expr[rownames(limma_out$coefficients), ge_data$pheno %>% filter(!!rlang::sym(contrast_info$variable) == contrast_info$active) %>% pull(sample_name)] %>% t %>% as.data.frame %>%
      lapply(function (x) {
        x[!is.na(x)] %>% length
      }) %>% unlist %>% unname

    diff_data_out$reference_n <- ge_data$expr[rownames(limma_out$coefficients), ge_data$pheno %>% filter(!!rlang::sym(contrast_info$variable) == contrast_info$reference) %>% pull(sample_name)] %>% t %>% as.data.frame %>%
      lapply(function (x) {
        x[!is.na(x)] %>% length
      }) %>% unlist %>% unname
  }


  ## I'll need to package these items a little differently for the covariates:
  diff_data_all_vars <- list(
    gene = rownames(limma_out$coefficients),
    coeff = limma_out$coefficients,
    p_value = limma_out$p.value,
    cohen_d = limma_out$coefficients / sqrt(limma_out$s2.post), # As recommended by Gordon Smyth: https://support.bioconductor.org/p/71747/
    variance = limma_out$s2.post * limma_out$stdev.unscaled^2 # --||-- in https://support.bioconductor.org/p/70175/
  )


  ge_results <- structure(
    list(
      diff_data = diff_data_out,
      diff_data_all_vars = diff_data_all_vars,
      dataset_name = ge_object$dataset_name,
      technology = tech,
      design_matrix = mm
    ), class = "diff_results"
  )

  # })

  return(ge_results)
}



run_dgsva <- function (
  ge_object,
  gene_sets = NULL#,
  # translate_to_entrezid = NULL
) {

  #' `run_dgsva`, as in "run differential GSVA".
  #' Here using gsva and limma.
  #'
  #' *Could be passing in a concrete contrast argument instead with a more
  #' complex contrast once needed (when/if implementing, remember that these
  #' will need to be supplied with `make.names`'d variable names) - for now
  #' working with simple contrasts, so don't need to implement this.
  #' *This expression should involve all levels from the "active_levels" and
  #' "reference_levels" lists.
  #'
  #' Remember to pre-merge any synonymous levels before using OMA!


  if (!is(ge_object, "dge_ready_data")) { stop("`ge_object` needs to be of `dge_ready_data` class!") }
  # lapply(ge_object, check_ge_obj)


  gene_sets_for_getter <- ifelse(is.null(gene_sets), "NULL", paste0('c("', paste0(gene_sets, collapse = '", "'), '")'))
  oma_tempscript <- file.path(tempdir(), "OMA-tempscript.R")
  oma_tempdata <- file.path(tempdir(), "OMA-tempdata.rds")
  geneset_getter_run <- paste0('library(OMA)
  out <- geneset_get(', gene_sets_for_getter, ')
  saveRDS(out, "', oma_tempdata, '")
  ')
  writeLines(geneset_getter_run, oma_tempscript)
  system(
    paste0("Rscript --vanilla --max-ppsize=500000 ", oma_tempscript)
  )
  gene_sets <- readRDS(oma_tempdata)

  ge_info <- ge_object$info
  contrast_info <- ge_info$contrast

  if (ge_info$contrast_type == "binary") { ## this is the case only for binary contrasts
    active_nicename <- make.names(contrast_info[["active"]])
    reference_nicename <- make.names(contrast_info[["reference"]])
    contrast_nicename <- paste0(active_nicename, "-", reference_nicename)
  }


  # diff_results <- lapply(1:length(ge_object$data), function (d_num) {

  tech <- ge_info$technology

  diff_data <- ge_object$data


  ## Right up front, `model.matrix` is not able to automatically handle any rows with NA's.
  nonna_rows <- apply(diff_data$pheno[, c(contrast_info$variable, diff_data$covariates), drop = FALSE], 1, function(x) all(!is.na(x)))
  diff_data$expr <- diff_data$expr[, nonna_rows]
  diff_data$pheno <- diff_data$pheno[nonna_rows, ]


  # if (is.null(translate_to_entrezid)) {
  #   x_to_entrezid <- rownames(diff_data$expr)
  # } else {
  #   x_to_entrezid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(diff_data$expr), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  # }


  # diff_data$for_gsva <- tibble(entrezid = x_to_entrezid) %>% bind_cols(diff_data$expr %>% as_tibble) %>% tidyr::drop_na(entrezid) %>% tibble::column_to_rownames("entrezid")
  diff_data$for_gsva <- fix_missing_tx(diff_data$expr, rownames(diff_data$expr), colnames(diff_data$expr))
  # diff_data$for_gsva <- tibble(entrezid = x_to_entrezid) %>% 
  #   bind_cols(diff_data$expr %>% as_tibble) %>% 
  #   tidyr::drop_na(entrezid) %>% 
  #   tibble::column_to_rownames("entrezid") %>% 
  #   apply(1, function (x) {
  #     if (all(!is.na(x))) {
  #       return(x)  # Return the row as is if all values are NA
  #     } else if (all(is.na(x))) {
  #       return(rep(0, length(x)))  # Return NA if all values are NA
  #     }
  #     na.replace <- median(x, na.rm = TRUE)  # Calculate median excluding NAs
  #     replace(x, is.na(x), na.replace)  # Replace NA with median
  #   }) %>% t


  # diff_data$gsva <- GSVA::gsva(diff_data$for_gsva %>% as.matrix, gene_sets, min.sz = 10, max.sz = 500, method = "gsva", rnaseq = TRUE, mx.diff = TRUE, verbose = FALSE)
  if (tech == "RNA-seq") {
    diff_data$gsva <- GSVA::gsva(diff_data$for_gsva, gene_sets, min.sz = 10, max.sz = 500, method = "gsva", kcdf = "Poisson", mx.diff = TRUE, verbose = FALSE)
  } else {
    diff_data$gsva <- GSVA::gsva(diff_data$for_gsva, gene_sets, min.sz = 10, max.sz = 500, method = "gsva", mx.diff = TRUE, verbose = FALSE)
  }

  diff_f <- ifelse(
    ge_info$contrast_type == "binary",
    paste0('~ 0 + diff_data[["pheno"]][, "', contrast_info$variable, '", drop = TRUE]'),
    paste0('~ diff_data[["pheno"]][, "', contrast_info$variable, '", drop = TRUE]')
  )
  if (!is.null(diff_data$covariates))
    diff_f <- paste0(diff_f, " + ", paste0(sapply(diff_data$covariates, function (x) paste0('diff_data$pheno[, "', x, '", drop = TRUE]')), collapse = " + "))
  diff_f <- as.formula(diff_f)

  # print(diff_data[["pheno"]][, contrast_info$variable])
  # print(contrast_info)
  # print(diff_data[["pheno"]][, "gender"])

  mm <- model.matrix(diff_f)
  if (ge_info$contrast_type == "binary") ## this is the case only for binary contrasts
    colnames(mm)[1:length(levels(diff_data[["pheno"]][, contrast_info$variable, drop = TRUE]))] <- levels(diff_data[["pheno"]][, contrast_info$variable, drop = TRUE])
  colnames(mm) <- make.names(colnames(mm))

  if (isTRUE(ge_info$use_sva)) {
    n_sv <- sva::num.sv(diff_data$gsva, mm, method = ge_info$sva_method)
    if (n_sv == 0)
      stop("Can't run with automated SVA - number of surrogate variables picked by SVA is 0!")
    sva_obj <- sva::sva(diff_data$gsva, mm, n.sv = n_sv)
    sva_sv <- sva_obj$sv
    colnames(sva_sv) <- paste0("sv", seq(1, n_sv))
    mm <- cbind(mm, sva_sv)
  }

  print(colnames(mm))

  if (ge_info$contrast_type == "binary") ## this is the case only for binary contrasts
    con <- limma::makeContrasts(contrasts = contrast_nicename, levels = mm)


  lmfit_partial <- purrr::safely(limma::lmFit) %>% purrr::partial(diff_data$gsva, design = mm)
  lmfit_realized <- lmfit_partial()$result

  if (!is.null(diff_data$block)) {
    corfit <- limma::duplicateCorrelation(diff_data$gsva, mm, block = diff_data$pheno[, diff_data$block, drop = TRUE])
    lmfit_partial_cor <- lmfit_partial %>% purrr::partial(block = diff_data$pheno[, diff_data$block, drop = TRUE], correlation = corfit$result$consensus.correlation)
    lmfit_partial_cor_realized <- lmfit_partial_cor()
    if (is.null(lmfit_partial_cor_realized$error)) {
      lmfit_realized <- lmfit_partial_cor_realized$result
    }
  }

  if (ge_info$contrast_type == "binary") ## this is the case only for binary contrasts
    limma_realized <- limma::contrasts.fit(
      lmfit_realized,
      con
    )
  else if (ge_info$contrast_type == "numeric")
    limma_realized <- lmfit_realized

  limma_out <- limma::eBayes(limma_realized)

  # limma::topTable(limma_out, number = Inf) %>% tibble::rownames_to_column(var = "gene") %>% head

  if (ge_info$contrast_type == "binary")
    ix_contrast <- 1
  else if (ge_info$contrast_type == "numeric")
    ix_contrast <- 2

  diff_data_out <- diff_data <- tibble(
    pathway = rownames(limma_out$coefficients),
    coeff = limma_out$coefficients[, ix_contrast] %>% as.vector(),
    p_value = limma_out$p.value[, ix_contrast] %>% as.vector(),
    cohen_d = limma_out$coefficients[, ix_contrast] / sqrt(limma_out$s2.post),
    variance = limma_out$s2.post * limma_out$stdev.unscaled[, ix_contrast]^2 # According to Gordon Smyth in https://support.bioconductor.org/p/70175/
  )


  if (ge_info$contrast_type == "binary") { ## this is the case only for binary contrasts
    diff_data_out$active_n <- diff_data$gsva[rownames(limma_out$coefficients), diff_data$pheno %>% filter(!!rlang::sym(contrast_info$variable) == contrast_info$active) %>% pull(sample_name)] %>% t %>% as.data.frame %>%
      lapply(function (x) {
        x[!is.na(x)] %>% length
      }) %>% unlist %>% unname

    diff_data_out$reference_n <- diff_data$gsva[rownames(limma_out$coefficients), diff_data$pheno %>% filter(!!rlang::sym(contrast_info$variable) == contrast_info$reference) %>% pull(sample_name)] %>% t %>% as.data.frame %>%
      lapply(function (x) {
        x[!is.na(x)] %>% length
      }) %>% unlist %>% unname
  }


  ## I'll need to package these items a little differently for the covariates:
  diff_data_all_vars <- list(
    pathway = rownames(limma_out$coefficients),
    coeff = limma_out$coefficients,
    p_value = limma_out$p.value,
    cohen_d = limma_out$coefficients / sqrt(limma_out$s2.post), # As recommended by Gordon Smyth: https://support.bioconductor.org/p/71747/
    variance = limma_out$s2.post * limma_out$stdev.unscaled^2 # --||-- in https://support.bioconductor.org/p/70175/
  )


  diff_results <- structure(
    list(
      diff_data = diff_data_out,
      diff_data_all_vars = diff_data_all_vars,
      dataset_name = ge_object$dataset_name,
      technology = tech,
      design_matrix = mm
    ), class = "diff_results"
  )

  # })

  return(diff_results)
}



run_ma <- function (
  diff_results,
  id_var = "gene",
  es_var = "coeff",
  variance_var = NULL,
  se_var = NULL,
  p_adj_method = c("BH",
  p.adjust.methods),
  Q_p_cutoff = 0.05,
  inclusion_cutoff = 1,
  parallel = FALSE,
  n_cores = future::availableCores() - 1
) {

  #' `run_ma` uses metafor's `rma` and `rma.mv` for meta-analysis.
  #'
  #' @param `dat_results` must be a data frame with, at the minimum, a "gene"
  #'   column, an effect size column, and a variance column.
  #' @param `dat_levels` needs to be passed if we're dealing with multilevel
  #'   data. If this is the case, the level has to be specified for each data-
  #'   set. This will trigger the calculation of multi-level fixed and random
  #'   effects.
  #' @param `es_var` here is unrestricted, but if the results are from `run_dge`,
  #'   then this should be either "cohen_d", or "coeff"; "coeff" is favored here,
  #'   as using "cohen_d" provides rescaling that also rescales the variance,
  #'   which is not accounted for by this method.
  #' @param `variance_var` here is unrestricted, but if the results are from
  #'   `run_dge`, then this should be "p_pvalue".
  #' @param `se_var` supply the name of the column containing the standard
  #'   errors if you've standard errors at your disposal rather than variances.
  #' @param `parallel` - if set to TRUE, furrr will be used instead of purrr in
  #'   the M-A loop

  p_adj_method <- match.arg(p_adj_method)


  if (!is.list(diff_results)) { stop("`diff_results` needs to be a list of differential results objects!") }
  if (!all(sapply(diff_results, function (x) is(x, "diff_results"))))


  if (!is.null(variance_var) & !is.null(se_var)) { stop("One of `variance_var` or `se_var` should be set, but not both!") }
  if (is.null(variance_var) & is.null(se_var)) variance_var <- "variance"

  se_or_variance <- ifelse(is.null(se_var), "variance", "se")
  sevar <- ifelse(is.null(se_var), variance_var, se_var)


  if (!(parallel %in% c(TRUE, FALSE))) { stop("`parallel` can only be set to TRUE or FALSE.") }
  if ((n_cores %% 1) != 0) { stop("`n_cores` must be an integer number!") }

  if (parallel) {
    if (.Platform$OS.type == "windows") {
      future::plan(future::multisession, workers = n_cores)
    } else {
      future::plan(future::multicore, workers = n_cores)
    }
    map <- furrr::future_map
  } else {
    map <- lapply
  }


  dat_levels <- sapply(diff_results, function (x) x$technology)


  multilevel <- length(unique(dat_levels)) > 1


  ma_info <- list(
    run_date = date(),
    es_variable = es_var,
    se_or_variance = se_or_variance,
    variance_variable = variance_var,
    se_variable = se_var,
    p_adj_method = p_adj_method,
    Q_p_cutoff = Q_p_cutoff,
    inclusion_cutoff = inclusion_cutoff
  )


  dat_names <- sapply(diff_results, function (x) x$dataset_name)


  if (is.null(dat_levels)) {
    dat_levels_ <- rep(NA, length(dat_names))
  } else {
    dat_levels_ <- dat_levels
  }


  dat_levels_df <- tibble(dataset = dat_names) %>% mutate(dat_level = dat_levels_)


  dat_results <- lapply(diff_results, function (x) x$diff_data)
  names(dat_results) <- dat_names


  ids_across <- dat_results %>% lapply(function (x) x %>% rename(idvar = !!id_var) %>% pull(idvar)) %>% purrr::reduce(union)

  idwise_es <- c(tibble(idvar = ids_across) %>% list, dat_results %>% map(function (x) x %>% rename(idvar = !!id_var) %>% select(idvar, !!rlang::sym(es_var)))) %>% purrr::reduce(left_join, by = "idvar")
  colnames(idwise_es) <- c("idvar", names(dat_results))
  # idwise_es_n <- idwise_es %>% tidyr::gather("key", "es", -gene) %>% select(-key) %>% tidyr::drop_na() %>% group_by(gene) %>% summarize(n = n()) %>% ungroup

  idwise_sevar <- c(tibble(idvar = ids_across) %>% list, dat_results %>% map(function (x) x %>% rename(idvar = !!id_var) %>% select(idvar, !!rlang::sym(sevar)))) %>% purrr::reduce(left_join, by = "idvar")
  colnames(idwise_sevar) <- c("idvar", names(dat_results))
  # idwise_var_n <- idwise_var %>% tidyr::gather("key", "var", -gene) %>% select(-key) %>% tidyr::drop_na() %>% group_by(gene) %>% summarize(n = n()) %>% ungroup
 if (!all((idwise_es %>% pull(idvar)) == (idwise_sevar %>% pull(idvar))))
    stop("All id's and all numbers of datasets have to be equal in the ES and variance outputs!")


  ma_ready <- idwise_es %>% tidyr::gather("dataset", "es", -idvar) %>% left_join(idwise_sevar %>% tidyr::gather("dataset", "sevar", -idvar), by = c("idvar", "dataset")) %>% tidyr::drop_na() %>% left_join(dat_levels_df, by = "dataset") %>% group_split(idvar)
  names(ma_ready) <- sapply(ma_ready, function (x) x[1, "idvar"])


  print("length(ma_ready):")
  print(length(ma_ready))
  if (multilevel) {
    ma_ready <- ma_ready %>% lapply(function (ma_set) {
      if ((ma_set$dat_level %>% unique %>% length) < 2)
        return(NULL)
      return(ma_set)
    })
    ma_ready[sapply(ma_ready, is.null)] <- NULL
  }
  print("length(ma_ready):")
  print(length(ma_ready))
  # stop()


  multilevel_modeling <- function (partial_model, random_bit, data) {
    if ((data$dat_level %>% unique %>% length) < 2) { stop("Two or more `dat_level`s necessary - condition not met in data!") }
    partial_model(random = random_bit, dat = data)
  }


  # ma_ready <- ma_ready[1:2] # TESTING
  ma_start_time <- Sys.time()

  if (!multilevel) {
    ma_results <- map(ma_ready, function (ma_set) {
      list(
        ma_fixed = purrr::safely(metafor::rma)(yi = ma_set$es, vi = ifelse(se_or_variance == "variance", ma_set$sevar, sqrt(ma_set$sevar)), method = "FE"),
        ma_random = purrr::safely(metafor::rma)(yi = ma_set$es, vi = ifelse(se_or_variance == "variance", ma_set$sevar, sqrt(ma_set$sevar)), method = "DL")
      )
    })
  } else {
    ma_results <- map(ma_ready, function (ma_set) {
      list(
        # ma_fixed = purrr::safely(multilevel_modeling)(partial_model = purrr::partial(metafor::rma.mv, yi = ma_set$es, V = ifelse(se_or_variance == "variance", ma_set$sevar, sqrt(ma_set$sevar))), random_bit = ~ 1 | dat_level, data = ma_set),# method = "DL"),
        ma_fixed = purrr::safely(metafor::rma.mv)(yi = ma_set$es, V = ifelse(se_or_variance == "variance", ma_set$sevar, sqrt(ma_set$sevar)), random = ~ 1 | dat_level, dat = ma_set),# method = "DL"),
        # ma_random = purrr::safely(multilevel_modeling)(partial_model = purrr::partial(metafor::rma.mv, yi = ma_set$es, V = ifelse(se_or_variance == "variance", ma_set$sevar, sqrt(ma_set$sevar))), random_bit = ~ 1 | dat_level/dataset, data = ma_set)#, method = "DL")
        ma_random = purrr::safely(metafor::rma.mv)(yi = ma_set$es, V = ifelse(se_or_variance == "variance", ma_set$sevar, sqrt(ma_set$sevar)), random = ~ 1 | dat_level/dataset, dat = ma_set)#, method = "DL")
      )
    })
  }

  ma_done_time <- Sys.time()
  ma_time_taken <- ma_done_time - ma_start_time

  cat("\nPure meta-analysis took:\n")
  print(ma_time_taken)


  ma_object <- tibble(idvar = names(ma_ready)) %>% bind_cols(
    tibble(
      n_datasets = sapply(ma_ready, nrow),
      datasets = sapply(ma_ready, function (x) x$dataset %>% list)
    )
  ) %>% bind_cols(
    ma_results %>% lapply(function (x) {

      if (is.null(x$ma_fixed$error)) {
        fixed_bit <- tibble(fixed_es = x$ma_fixed$result$beta[1], fixed_se = x$ma_fixed$result$se, fixed_ci_lb = x$ma_fixed$result$ci.lb, fixed_ci_ub = x$ma_fixed$result$ci.ub, fixed_pval = x$ma_fixed$result$pval)
      } else {
        fixed_bit <- tibble(fixed_es = NA, fixed_se = NA, fixed_ci_lb = NA, fixed_ci_ub = NA, fixed_pval = NA)
      }

      if (is.null(x$ma_random$error)) {
        random_bit <- tibble(random_es = x$ma_random$result$beta[1], random_se = x$ma_random$result$se, random_ci_lb = x$ma_random$result$ci.lb, random_ci_ub = x$ma_random$result$ci.ub, random_pval = x$ma_random$result$pval, tau2 = x$ma_random$result$tau2, QEp = x$ma_random$result$QEp)
      } else {
        random_bit <- tibble(random_es = NA, random_se = NA, random_ci_lb = NA, random_ci_ub = NA, random_pval = NA, tau2 = NA, QEp = NA)
      }

      bind_cols(fixed_bit, random_bit)

    }) %>% bind_rows
  ) %>% filter(!(is.na(fixed_es) & is.na(random_es)), n_datasets >= inclusion_cutoff) %>% mutate(
    fixed_adj_pval = p.adjust(fixed_pval, method = p_adj_method),
    random_adj_pval = p.adjust(random_pval, method = p_adj_method),
    hybrid_adj_pval = ifelse(!is.na(QEp) & (QEp < Q_p_cutoff), random_pval, fixed_pval) %>% p.adjust(method = p_adj_method),
    hybrid_es = ifelse(!is.na(QEp) & (QEp < Q_p_cutoff), random_es, fixed_es),
    hybrid_se = ifelse(!is.na(QEp) & (QEp < Q_p_cutoff), random_se, fixed_se)
  ) ## ma_object

  filtered_adj_p <- inclusion_cutoff:max(ma_object$n_datasets) %>% lapply(function (cutoff) {
    ma_object %>% filter(n_datasets >= cutoff) %>% select(idvar, fixed_pval, random_pval, QEp) %>%
      transmute(
        idvar,
        fixed_adj_pval = p.adjust(fixed_pval, method = p_adj_method),
        random_adj_pval = p.adjust(random_pval, method = p_adj_method),
        hybrid_adj_pval = ifelse(!is.na(QEp) & (QEp <= Q_p_cutoff), random_pval, fixed_pval) %>% p.adjust(method = p_adj_method)
      )
  })
  filtered_adj_p <- list("fixed_adj_pval", "random_adj_pval", "hybrid_adj_pval") %>% lapply(function (pval) {
    filtered_adj_p %>% purrr::reduce(function (x, y) {left_join(x, y %>% select(!!id_var := idvar, !!rlang::sym(pval)), by = id_var)}, .init = ma_object %>% select(!!id_var := idvar))
  })
  names(filtered_adj_p) <- c("fixed_adj_pval", "random_adj_pval", "hybrid_adj_pval")

  return(
    list(
      ma = ma_object %>% rename(!!id_var := idvar),
      info = ma_info,
      adj_pval = filtered_adj_p,
      adj_pval_cutoffs = inclusion_cutoff:max(ma_object$n_datasets)
    )
  )

}
