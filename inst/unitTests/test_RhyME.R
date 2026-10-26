test_RhyME <- function() {
    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    fit <- fitCosinor(example$se[1:5, ], subject_col = "subject",
        prior = "none", BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    result <- rowData(fit)[, startsWith(names(rowData(fit)), "RhyME_")]

    checkEquals(dim(result), c(5L, 13L))
    checkEquals(names(result), paste0("RhyME_", c(
        "b0", "b_sin", "b_cos", "log_phi", "tau", "score_chisq", "pvalue",
        "padj", "logit_amplitude", "amplitude", "phase", "mesor", "converged")))
    checkTrue(all(is.finite(result$RhyME_b0)))
    checkTrue(all(is.finite(result$RhyME_pvalue)))
    checkTrue(all(result$RhyME_pvalue >= 0 & result$RhyME_pvalue <= 1))
    checkTrue(all(result$RhyME_converged))
    checkTrue(inherits(plotCosinorFit(fit), "ggplot"))
}

test_registered_backend <- function() {
    previous <- BiocParallel::registered()
    on.exit({
        for (param in rev(previous)) BiocParallel::register(param)
    }, add = TRUE)
    withr::local_options(list(SerialParam = NULL))
    log_dir <- tempfile("RhyME-bpparam-")
    dir.create(log_dir)
    on.exit(unlink(log_dir, recursive = TRUE), add = TRUE)
    BiocParallel::register(BiocParallel::SerialParam(
        log = TRUE, logdir = log_dir))

    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    se <- example$se[1:5, ]
    fit_default <- fitCosinor(se, subject_col = "subject",
        prior = "none", verbose = FALSE)
    logs <- list.files(log_dir, full.names = TRUE)
    checkTrue(length(logs) > 0L)

    unlink(logs)
    fit_explicit <- fitCosinor(se, subject_col = "subject",
        prior = "none", BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    checkEquals(length(list.files(log_dir)), 0L)
    checkEquals(as.data.frame(rowData(fit_default)),
        as.data.frame(rowData(fit_explicit)))
}

test_row_annotations_preserved <- function() {
    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    se <- example$se[1:5, ]
    reference <- fitCosinor(se, subject_col = "subject", prior = "none",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    result_names <- setdiff(names(rowData(reference)), names(rowData(se)))

    # Include all generic names previously removed by the result assembler.
    legacy <- c(sub("^RhyME_", "", result_names), "RhyME",
        "status_null", "status_full", "error_null", "error_full", "phi",
        "amplitude_link", "amplitude_response", "mesor_link", "mesor_response")
    annotations <- S4Vectors::DataFrame(lapply(
        stats::setNames(legacy, legacy), function(x) rep(NA_real_, nrow(se))))
    annotations$converged <- rep(FALSE, nrow(se))
    annotations$gene <- factor(letters[seq_len(nrow(se))])
    annotations$RhyME_note <- rep("Keep this annotation", nrow(se))
    rownames(annotations) <- rownames(se)
    rowData(se) <- annotations

    fit <- fitCosinor(se, subject_col = "subject", prior = "none",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    checkEquals(as.data.frame(rowData(fit)[, names(annotations), drop = FALSE]),
        as.data.frame(annotations))
    checkEquals(as.data.frame(rowData(fit)[, result_names]),
        as.data.frame(rowData(reference)[, result_names]))
    checkEquals(S4Vectors::metadata(fit)$RhyME$fitSummary,
        S4Vectors::metadata(reference)$RhyME$fitSummary)

    # Missing generic p-values and coefficients must not affect the plot.
    checkEquals(ggplot2::ggplot_build(plotCosinorFit(fit))$data,
        ggplot2::ggplot_build(plotCosinorFit(reference))$data)
    selected <- rev(rownames(fit)[1:2])
    checkEquals(ggplot2::ggplot_build(plotCosinorFit(fit, cpgs = selected))$data,
        ggplot2::ggplot_build(plotCosinorFit(reference, cpgs = selected))$data)
}

test_refit_replaces_only_results <- function() {
    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    se <- example$se[1:5, ]
    rowData(se)$pvalue <- seq_len(nrow(se)) / nrow(se)
    rowData(se)$RhyME_note <- letters[seq_len(nrow(se))]
    fit <- fitCosinor(se, subject_col = "subject", prior = "none",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    expected <- rowData(fit)
    result_names <- setdiff(names(expected), names(rowData(se)))
    for (name in result_names) {
        rowData(fit)[[name]] <- rep(
            if (name == "RhyME_converged") FALSE else -1, nrow(fit))
    }

    refit <- fitCosinor(fit, subject_col = "subject", prior = "none",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    checkEquals(as.data.frame(rowData(refit)), as.data.frame(expected))
    checkEquals(anyDuplicated(names(rowData(refit))), 0L)
}

test_plot_requires_prefixed_results <- function() {
    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    fit <- fitCosinor(example$se[1:5, ], subject_col = "subject",
        prior = "none", BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    names(rowData(fit)) <- sub("^RhyME_", "", names(rowData(fit)))
    error <- tryCatch(plotCosinorFit(fit), error = conditionMessage)
    checkTrue(is.character(error))
    checkTrue(grepl("Refit with fitCosinor().", error, fixed = TRUE))
}
