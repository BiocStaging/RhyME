#' Plot Fitted Cosinor Curves
#'
#' Plots observed beta values and fitted marginal mean curves over the
#' observed time range, using the design saved by `fitCosinor()`. Curves
#' describe the mean covariate profile; points are not covariate-adjusted.
#'
#' @param se A `SummarizedExperiment` returned by `fitCosinor()`.
#' @param cpgs CpG identifiers in panel order. `NULL` selects up to six CpGs
#' with usable estimates, ordered by p-value.
#' @param ncol Number of panels per row.
#'
#' @return A `ggplot2` object.
#' @seealso [fitCosinor()]
#' @export
#'
#' @examples
#' ex <- readRDS(system.file("extdata", "example_data.rds", package = "RhyME"))
#' fit <- fitCosinor(ex$se[1, ], subject_col = "subject", prior = "none",
#'     verbose = FALSE)
#' plotCosinorFit(fit)
plotCosinorFit <- function(se, cpgs = NULL, ncol = 3) {
    res <- rowData(se)
    design <- metadata(se)$RhyME$design
    if (is.null(design)) {
        stop("'se' does not contain RhyME fit metadata.", call. = FALSE)
    }
    period <- design$period

    if (is.null(cpgs)) {
        eligible <- which(!is.na(res$pvalue) & !is.na(res$b0))
        eligible <- eligible[order(res$pvalue[eligible])]
        cpgs <- rownames(res)[utils::head(eligible, 6L)]
    }

    if (!length(cpgs) || anyNA(cpgs) || any(!cpgs %in% rownames(res)))
        stop("Select existing CpG names with usable parameter estimates.", call. = FALSE)
    if (anyNA(res[cpgs, "b0", drop = TRUE]))
        stop("A requested CpG has no usable full-model parameter estimates.", call. = FALSE)
    input <- .getRhyMEInput(se[cpgs, ], design$assay_name, design$time_col)
    times <- input$col_data[[design$time_col]]
    beta_mat <- input$beta

    raw_list <- lapply(cpgs, function(cg) {
        data.frame(
            cpg  = cg,
            time = times,
            beta = as.numeric(beta_mat[cg, ]),
            stringsAsFactors = FALSE
        )
    })
    raw_df <- do.call(rbind, raw_list)

    tgrid <- sort(unique(c(times, seq(min(times), max(times), length.out = 200))))
    omega <- 2 * pi / period
    gh <- .normalGH(15L)
    mean_curve <- function(eta, tau) {
        if (identical(design$mode, "mixed") && is.finite(tau) && tau > 0)
            drop(crossprod(gh$weights,
                stats::plogis(outer(tau * gh$nodes, eta, "+"))))
        else stats::plogis(eta)
    }
    fit_list <- lapply(cpgs, function(cg) {
        b0 <- res[cg, "b0", drop = TRUE]
        bs <- res[cg, "b_sin", drop = TRUE]
        bc <- res[cg, "b_cos", drop = TRUE]
        eta <- b0 + bs * sin(omega * tgrid) + bc * cos(omega * tgrid)
        tau <- res[cg, "tau", drop = TRUE]
        data.frame(
            cpg  = cg,
            time = tgrid,
            fitted = mean_curve(eta, tau),
            stringsAsFactors = FALSE
        )
    })
    fit_df <- do.call(rbind, fit_list)
    cpg_order <- unique(cpgs)
    raw_df$cpg <- factor(raw_df$cpg, levels = cpg_order)
    fit_df$cpg <- factor(fit_df$cpg, levels = cpg_order)

    ggplot2::ggplot() +
        ggplot2::geom_point(
            data = raw_df,
            ggplot2::aes(x = .data$time, y = .data$beta),
            alpha = 0.4, size = 1
        ) +
        ggplot2::geom_line(
            data = fit_df,
            ggplot2::aes(x = .data$time, y = .data$fitted),
            color = "red", linewidth = 0.8
        ) +
        ggplot2::facet_wrap(~ cpg, ncol = ncol, scales = "free_y") +
        ggplot2::labs(x = "Time (hours)", y = "Beta value",
                      title = "Cosinor fit") +
        ggplot2::theme_bw()
}
