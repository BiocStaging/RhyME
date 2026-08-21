#' Plot fitted cosinor curves for selected CpGs
#'
#' Overlays the fitted cosinor curve on raw beta values for a set of CpGs.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] object
#'   returned by [fitCosinor()].
#' @param cpgs Character vector of CpG identifiers to plot. If NULL, the top
#'   6 CpGs by p-value are plotted.
#' @param assay_name Character or NULL, assay name in `se`. The default NULL
#'   uses the first assay in `se`.
#' @param time_col Character, time column in colData (default `"time"`).
#' @param ncol Integer, number of columns in facet layout (default 3).
#' @param n_grid Integer, resolution of fitted curve grid (default 200).
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' ex <- readRDS(system.file("extdata", "example_data.rds",
#'                          package = "RhyME"))
#' se <- ex$se[seq_len(5), ]
#' fit <- fitCosinor(se, time_col = "time", subject_col = "subject",
#'                   prior = "none", verbose = FALSE)
#' plotCosinorFit(fit, cpgs = rownames(fit)[1], ncol = 1)
plotCosinorFit <- function(se, cpgs = NULL,
                           assay_name = NULL,
                           time_col = "time",
                           ncol = 3, n_grid = 200) {

    res <- rowData(se)
    design <- metadata(se)$RhyME$design
    period <- design$period

    if (is.null(cpgs)) {
        top <- topCpGs(se, n = 6)
        cpgs <- top$cpg
    }

    input <- .getRhyMEInput(se, assay_name, time_col)
    times <- input$col_data[[time_col]]
    beta_mat <- input$beta

    # Build raw data
    raw_list <- lapply(cpgs, function(cg) {
        data.frame(
            cpg  = cg,
            time = times,
            beta = as.numeric(beta_mat[cg, ]),
            stringsAsFactors = FALSE
        )
    })
    raw_df <- do.call(rbind, raw_list)

    # Build fitted curves
    tgrid <- seq(0, period, length.out = n_grid)
    omega <- 2 * pi / period
    fit_list <- lapply(cpgs, function(cg) {
        b0 <- res[cg, "b0", drop = TRUE]
        bs <- res[cg, "b_sin", drop = TRUE]
        bc <- res[cg, "b_cos", drop = TRUE]
        if (is.na(b0)) return(NULL)
        eta <- b0 + bs * sin(omega * tgrid) + bc * cos(omega * tgrid)
        data.frame(
            cpg  = cg,
            time = tgrid,
            fitted = .invlogit(eta),
            stringsAsFactors = FALSE
        )
    })
    fit_df <- do.call(rbind, fit_list)

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
