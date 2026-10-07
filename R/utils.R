.getRhyMEInput <- function(se, assay_name, time_col, subject_col = NULL,
                           covariates = NULL) {
    if (!inherits(se, "SummarizedExperiment"))
        stop("'se' must be a SummarizedExperiment.", call. = FALSE)
    if (!nrow(se) || !ncol(se) || is.null(rownames(se)) ||
        anyNA(rownames(se)) || any(!nzchar(rownames(se))) ||
        anyDuplicated(rownames(se)))
        stop("Input needs at least one CpG and sample, with unique CpG names.",
             call. = FALSE)
    for (name in list(time_col, subject_col, assay_name)) {
        if (!is.null(name) && (!is.character(name) || length(name) != 1L ||
                              is.na(name) || !nzchar(name)))
            stop("Column and assay names must be single nonempty strings.",
                 call. = FALSE)
    }
    if (is.null(time_col)) stop("'time_col' is required.", call. = FALSE)
    cd <- as.data.frame(colData(se))
    missing <- setdiff(c(time_col, subject_col, covariates), names(cd))
    if (length(missing))
        stop("Missing colData columns: ", paste(missing, collapse = ", "),
             call. = FALSE)
    if (!is.numeric(cd[[time_col]]))
        stop("Collection times must be numeric hours.", call. = FALSE)
    if (!is.null(subject_col)) {
        ids <- cd[[subject_col]]
        if (anyNA(ids) || any(!nzchar(as.character(ids))) ||
            (is.numeric(ids) && any(!is.finite(ids))))
            stop("Subject identifiers must be complete and finite.", call. = FALSE)
    }
    if (is.null(assay_name)) {
        beta <- assay(se, 1L)
        available <- assayNames(se)
        if (length(available) && !is.na(available[1]) && nzchar(available[1]))
            assay_name <- available[1]
    } else {
        beta <- assay(se, assay_name)
    }
    if (!is.matrix(beta) || !is.numeric(beta))
        stop("The assay must be a numeric matrix.", call. = FALSE)
    if (any(beta < 0 | beta > 1, na.rm = TRUE))
        stop("Observed beta values must be within [0, 1]; missing values may be NA.",
             call. = FALSE)
    keep <- is.finite(cd[[time_col]])
    list(beta = beta[, keep, drop = FALSE], col_data = cd[keep, , drop = FALSE],
         assay_name = assay_name)
}

.makeCovariateMatrix <- function(col_data, covariates = NULL) {
    if (!length(covariates))
        return(matrix(numeric(), nrow = nrow(col_data), ncol = 0L))
    if (!is.character(covariates) || anyNA(covariates) ||
        any(!nzchar(covariates)) || anyDuplicated(covariates))
        stop("'covariates' must contain unique colData column names.", call. = FALSE)
    cov_data <- col_data[, covariates, drop = FALSE]
    if (any(!stats::complete.cases(cov_data)))
        stop("Adjustment covariates cannot contain missing values.", call. = FALSE)
    X <- stats::model.matrix(~ ., data = cov_data)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    X <- sweep(X, 2L, colMeans(X), FUN = "-")
    colnames(X) <- paste0("cov_", make.names(colnames(X), unique = TRUE))
    X
}

.checkPQLDesign <- function(X, subject = NULL) {
    if (nrow(X) <= ncol(X) || qr(X)$rank != ncol(X))
        stop("A full-rank design and more observations than fixed-effect coefficients are required.")
    if (!is.null(subject)) {
        counts <- table(subject)
        if (length(counts) < 2L || !any(counts > 1L))
            stop("Mixed PQL requires multiple subjects with repeated observations.")
    }
    invisible(NULL)
}

.validatePQLPrior <- function(prior) {
    if (!is.numeric(prior$mean) || !is.numeric(prior$sd) ||
        length(prior$mean) != 1L || length(prior$sd) != 1L ||
        !is.finite(prior$mean) || !is.finite(prior$sd) || prior$sd <= 0)
        stop("A supplied prior must have a finite mean and positive finite sd on log(phi).",
             call. = FALSE)
    invisible(NULL)
}

.pqlBatches <- function(beta, workers) {
    if (workers == 1L) return(list(beta))
    index <- seq_len(nrow(beta))
    groups <- split(index, ceiling(index / ceiling(nrow(beta) / workers)))
    unname(lapply(groups, function(i) beta[i, , drop = FALSE]))
}

.pqlPilotIndex <- function(n) {
    if (n <= 10000L) return(seq_len(n))
    withr::with_seed(42, sort(sample.int(n, 10000L)))
}

.estimateLogPhiPrior <- function(log_phi) {
    valid <- log_phi[is.finite(log_phi)]
    if (length(valid) < 1000L)
        stop("EB prior estimation requires at least 1,000 usable null pilot fits.",
             call. = FALSE)
    list(mean = unname(stats::median(valid)),
         sd = max(unname(stats::IQR(valid)) / 1.349, 1e-3), n = length(valid))
}
