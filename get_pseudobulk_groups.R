

# helper to create pseudobulk from single-cell data
pseudobulk <- function(counts, grouping) {
  group_levels <- levels(grouping)
  has_group <- !is.na(grouping)
  if (inherits(x = counts, what = 'dgCMatrix')) {
    if (requireNamespace("sparseMatrixStats", quietly = TRUE)) {
      mat <- sapply(group_levels, function(gr) {
        sparseMatrixStats::rowSums2(x = counts, cols = has_group & grouping == gr)
      })
    } else {
      message('Consider installing the "sparseMatrixStats" package for faster pseudobulk aggregation')
      mat <- sapply(group_levels, function(gr) {
        Matrix::rowSums(x = counts[, has_group & grouping == gr, drop = FALSE])
      })
    }
  } else {
    mat <- t(rowsum(x = t(counts[, has_group, drop = FALSE]),
                    group = grouping[has_group]))
  }

  colnames(mat) <- group_levels
  rownames(mat) <- rownames(counts)
  return(mat)
}


make_pseudoreplicate_label <- function(grouping, G = 3) {
  tab <- table(grouping)
  if (any(tab < G)) {
    stop('too few cells in at least one group')
  }
  # create replicate labels per group
  replicate <- numeric(length = length(grouping))
  group_idx <- 1
  for (gl in levels(grouping)) {
    sel <- which(grouping == gl)
    replicate[sel] <- sample(1:tab[gl] %% G) + group_idx
    group_idx <- group_idx + G
  }

  # save replicate as csv for reproducibility
  write.csv(replicate, file = 'replicate4.csv')

  return(replicate)
}

# make pseudo-bulk from single-cell count data
make_pb <- function(counts, grouping, replicate_label) {
  grouping <- droplevels(as.factor(grouping))
  group_levels <- levels(grouping)

  if (is.null(replicate_label) | length(replicate_label) == 0) {
    replicate_label <- make_pseudoreplicate_label(grouping, G = 3)
  }
  replicate_label <- droplevels(as.factor(replicate_label))

  group_replicate <- interaction(grouping, replicate_label, drop = TRUE)

  # create pseudobulk
  pb <- pseudobulk(counts = counts, grouping = group_replicate)

  # save pseudobulk as matrix for reproducibility
  write.csv(pb, file = paste0('pseudobulk_cluster_', clust, '.csv'), row.names = TRUE, col.names = TRUE)

  # create pseudobulk meta data
  pb_md <- as.data.frame(table(grouping, replicate_label))
  pb_md <- pb_md[pb_md$Freq > 0, ]

  return(list(counts = pb, md = pb_md))
}

run_deseq_simple <- function(counts, grouping, order_results, lfc_shrinkage) {
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, data.frame(grouping = grouping), ~ grouping)
  dds <- DESeq2::DESeq(dds, test = 'Wald', quiet = TRUE)
  if (is.null(lfc_shrinkage)) {
    res <- DESeq2::results(dds)
  }
  else {
    if (lfc_shrinkage %in% c("apeglm", "ashr", "normal")) {
      coef_name <- DESeq2::resultsNames(dds)[2]
      res <- DESeq2::lfcShrink(dds, coef = coef_name, type = lfc_shrinkage)
    }
    else {
      stop('lfc_shrinkage should be set to NULL or one of "apeglm", "ashr", "normal"')
    }
  }
  res <- as.data.frame(res) %>%
    tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename('ave_expr' = 'baseMean', 'log_fc' = 'log2FoldChange')
  # some shrinkage methods drop the stat column - add NA
  if (!'stats' %in% colnames(res)) {
    res <- tibble::add_column(res, stat = NA_real_, .before = 'pvalue')
  }
  # select the columns we want to return
  res <- dplyr::select(res, .data$feature, .data$ave_expr, .data$log_fc, .data$stat, .data$pvalue, .data$padj)
  if (order_results) {
    res <- dplyr::arrange(res, .data$pvalue, -abs(.data$stat))
  }
  return(res)
}


run_de_simple <- function(counts, grouping, replicate_label, method, order_results, lfc_shrinkage, verbosity) {
  if (!any(method %in% c('deseq')) | length(method) != 1) {
    stop("method needs to be of length one and in c('deseq')")
  }

  # we are only dealing with bulk methods, so create pseudobulk data
  pb <- make_pb(counts, grouping, replicate_label)
  if (verbosity > 1) {
    message(paste0(c('group vs replicate table', utils::capture.output(pb$md)), collapse = "\n"))
  }

  if (verbosity > 0 && !is.null(lfc_shrinkage) && method == 'deseq') {
    message(sprintf('calculate shrunken log2 fold changes of type %s', lfc_shrinkage))
  }

  # the DE methods we use, treat the first group level as reference
  # we want this to be group2, so relevel here
  pb$md$grouping <- stats::relevel(pb$md$grouping, ref = levels(pb$md$grouping)[2])

  if (method == 'deseq') {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop('DESeq2 package not found - please install it to use "method = deseq"')
    }
    res <- run_deseq_simple(pb$counts, pb$md$grouping, order_results, lfc_shrinkage)
  }

  # add some feature-level stats
  det_rate <- detection_rate(counts, grouping)[res$feature, ]
  colnames(det_rate) <- c('rate1', 'rate2')
  rownames(det_rate) <- NULL
  res <- cbind(res, det_rate)

  return(res)
}


#' @importFrom magrittr %>%
#' @importFrom rlang .data
run_de_one_comp <- function(counts, grouping, replicate_label, comp, method, order_results, lfc_shrinkage, verbosity) {
  this_grouping <- as.numeric(grouping %in% comp[[3]]) + as.numeric(grouping %in% comp[[4]])*2
  this_grouping <- factor(this_grouping, levels = c(1, 2))
  levels(this_grouping) <- c(comp[[1]], comp[[2]])
  res <- run_de_simple(counts, this_grouping, replicate_label, method, order_results, lfc_shrinkage, verbosity) %>%
    dplyr::mutate(group1 = factor(comp[[1]]), group2 = factor(comp[[2]]))
  return(res)
}

#' @importFrom magrittr %>%
#' @importFrom rlang .data
run_de_comparisons <- function(counts, grouping, replicate_label, comparisons, method, order_results, lfc_shrinkage, verbosity) {
  # keep track of progress if progressr package is installed
  p <- function(x) {
    if (verbosity > 0) {
      message(x)
    }
    invisible()
  }
  if (requireNamespace("progressr", quietly = TRUE)) {
    p <- progressr::progressor(along = comparisons)
  }

  if (requireNamespace("future.apply", quietly = TRUE)) {
    res_list <- future.apply::future_lapply(comparisons, function(comp) {
      p(sprintf('Comparing %s vs %s', comp[[1]], comp[[2]]))
      run_de_one_comp(counts, grouping, replicate_label, comp, method, order_results, lfc_shrinkage, verbosity)
    }, future.seed = TRUE)
  } else {
    message('Consider installing the "future.apply" package for parallelization')
    res_list <- lapply(comparisons, function(comp) {
      p(sprintf('Comparing %s vs %s', comp[[1]], comp[[2]]))
      run_de_one_comp(counts, grouping, replicate_label, comp, method, order_results, lfc_shrinkage, verbosity)
    })
  }

  res <- dplyr::bind_rows(res_list)
  return(as.data.frame(res))
}

# helper to set up list of comparisons
# each comparison is a list
# name1, name2, labels grp1, labels grp2
set_up_comparisons <- function(group_levels, compare, compare_is_ref, verbosity) {
  n <- length(group_levels)
  if (n < 2) {
    stop('Did not find two or more group levels - make sure grouping variable is a factor with at least two levels')
  }

  if (compare[1] == 'each_vs_rest' && n == 2) {
    compare <- group_levels
    if (verbosity > 0) {
      message('There are only two groups in the data. Changing compare argument from "each_vs_rest" to group levels')
    }
  }

  if (compare[1] == 'each_vs_rest') {
    comparisons <- lapply(group_levels, function(x) list(x, paste0('not ', x), x, setdiff(group_levels, x)))
    return(comparisons)
  }

  if (compare[1] == 'all_vs_all') {
    comparisons <- list()
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        comparisons[[length(comparisons) + 1]] <- list(group_levels[i], group_levels[j], group_levels[i], group_levels[j])
      }
    }
    return(comparisons)
  }

  if (inherits(x = compare, what = 'character')) {
    if (length(compare) == 1) {
      if (!all(compare %in% group_levels)) {
        stop('Group 1 not found - please check your compare argument')
      }
      if (compare_is_ref) {
        ref <- compare[1]
        comparisons <- lapply(setdiff(group_levels, ref), function(x) list(x, ref, x, ref))
      } else {
        x <- compare[1]
        comparisons <- list(list(x, paste0('not ', x), x, setdiff(group_levels, x)))
      }
      return(comparisons)
    }
    if (length(compare) == 2) {
      if (!all(compare %in% group_levels)) {
        stop('Group 1 or Group 2 not found - please check your compare argument')
      }
      if (compare[1] == compare[2]) {
        stop('Group 1 and 2 need to be different - please check your compare argument')
      }
      comparisons <- list(list(compare[1], compare[2], compare[1], compare[2]))
      return(comparisons)
    }
    if (length(compare) > 2) {
      stop('Too many elements in compare argument')
    }
  }

  if (inherits(x = compare, what = 'list') &&
      length(compare) == 2 &&
      all(unlist(lapply(compare, inherits, what = 'character')))) {
    compare <- lapply(compare, unique)
    if (length(intersect(compare[[1]], compare[[2]])) > 0) {
      stop('Intersection between group 1 and 2 - please check your compare argument')
    }
    comparisons <- list(list('group1', 'group2', compare[[1]], compare[[2]]))
    return(comparisons)
  }

  stop("Make sure the compare argument is 'each_vs_rest' or 'all_vs_all' or a length 1 or 2
        character vector with entries present in the group_labels argument or
        a list of length 2 with each entry being a character vector of group labels")
}

print_comparisons <- function(comparisons, verbosity) {
  if (verbosity > 0) {
    message(sprintf('Doing %d comparisons', length(comparisons)))
  }
  if (verbosity > 1) {
    for (i in 1:length(comparisons)) {
      message(sprintf('%d:\t%s vs %s', i, comparisons[[i]][[1]], comparisons[[i]][[2]]))
    }
  }
  return(invisible())
}

# function to extract the count matrix and grouping factor from input
get_data <- function(object, meta_data, group_column, replicate_column, verbosity) {
  replicate_label <- NULL

  if (inherits(x = object, what = 'Seurat')) {
    if (!('RNA' %in% SeuratObject::Assays(object))) {
      stop('Did not find "RNA" assay in Seurat object')
    }
    if (!is.null(meta_data) & verbosity > 0) {
      message('input is Seurat object - the meta_data argument will be ignored')
    }
    counts <- SeuratObject::GetAssayData(object[["RNA"]], slot = "counts")
    meta_data <- object[[]]
    if (is.null(group_column)) {
      if (verbosity > 0) {
        message('group_column argument is not set - will use SeuratObject::Idents()')
      }
      grouping <- SeuratObject::Idents(object)
    } else {
      grouping <- meta_data[, group_column]
    }
  }

  if (inherits(x = object, what = 'SingleCellExperiment')) {
    if (!is.null(meta_data) & verbosity > 0) {
      message('input is SingleCellExperiment object - the meta_data argument will be ignored')
    }
    counts <- SingleCellExperiment::counts(object)
    meta_data <- as.data.frame(SummarizedExperiment::colData(object))
    if (is.null(group_column)) {
      if (verbosity > 0) {
        message('group_column argument is not set - will use SingleCellExperiment::colLabels()')
      }
      grouping <- SingleCellExperiment::colLabels(object)
    } else {
      grouping <- meta_data[, group_column]
    }
  }

  if (inherits(x = object, what = c('dgCMatrix', 'matrix'))) {
    counts <- object
    grouping <- meta_data[, group_column]
  }

  if (ncol(counts) != length(grouping)) {
    stop('Number of cells (columns) in count matrix does not match length of grouping
         variable - check meta_data and group_column arguments')
  }

  if (!is.null(replicate_column)) {
    replicate_label <- meta_data[, replicate_column]
  }

  grouping <- droplevels(as.factor(grouping))
  replicate_label <- droplevels(as.factor(replicate_label))

  if ((is.null(replicate_label) | length(replicate_label) == 0) & verbosity > 0) {
    message('No replicates specified - will assign cells randomly to pseudoreplicates')
  }

  return(list(counts = counts, grouping = grouping, replicate_label = replicate_label))
}


# for each feature determine fraction of non-zero cells per group
detection_rate <- function(counts, grouping) {
  group_levels <- levels(grouping)
  has_group <- !is.na(grouping)
  if (inherits(x = counts, what = 'dgCMatrix')) {
    if (requireNamespace("sparseMatrixStats", quietly = TRUE)) {
      mat <- sapply(group_levels, function(gr) {
        sparseMatrixStats::rowMeans2(x = counts > 0, cols = has_group & grouping == gr)
      })
    } else {
      message('Consider installing the "sparseMatrixStats" package for faster pseudobulk aggregation')
      mat <- sapply(group_levels, function(gr) {
        Matrix::rowMeans(x = counts[, has_group & grouping == gr, drop = FALSE] > 0)
      })
    }
  } else {
    tab <- table(grouping)
    mat <- t(rowsum(x = t((counts[, has_group, drop = FALSE] > 0) + 0),
                    group = grouping[has_group]) / as.numeric(tab))
  }

  colnames(mat) <- group_levels
  rownames(mat) <- rownames(counts)
  return(mat)
}


findDE <- function(object,
                   meta_data = NULL,
                   group_column = NULL,
                   replicate_column = NULL,
                   compare = 'each_vs_rest',
                   compare_is_ref = FALSE,
                   method = 'edger',
                   order_results = TRUE,
                   lfc_shrinkage = NULL,
                   verbosity = 1) {
  # extract the data from the input object
  de_data <- get_data(object, meta_data, group_column, replicate_column, verbosity)
  # set up the comparisons
  group_levels = levels(x = de_data$grouping)
  comparisons <- set_up_comparisons(group_levels = group_levels, compare = compare,
                                    compare_is_ref = compare_is_ref, verbosity = verbosity)
  print_comparisons(comparisons, verbosity)
  # run DE
  run_de_comparisons(counts = de_data$counts,
                     grouping = de_data$grouping,
                     replicate_label = de_data$replicate_label,
                     comparisons = comparisons,
                     method = method,
                     order_results = order_results,
                     lfc_shrinkage = lfc_shrinkage,
                     verbosity = verbosity)
}
