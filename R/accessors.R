#' Compute circadian parameters from cosinor coefficients
#'
#' @param b0 Numeric vector of intercepts.
#' @param b_sin Numeric vector of sin coefficients.
#' @param b_cos Numeric vector of cos coefficients.
#' @param period Numeric, the period in hours.
#' @param tau Numeric vector of random-intercept SD estimates. Used to
#'   recover marginal response-scale curves in mixed mode.
#' @param grid_by Numeric, step size for the numerical time grid.
#' @param mc_n Integer, number of fixed normal quantile nodes used to average
#'   over the random intercept distribution.
#' @return A list with amplitude_link, amplitude_response, phase,
#'   mesor_link, mesor_response.
#' @noRd
.computeCircadianParams <- function(b0, b_sin, b_cos, period = 24,
                                    tau = NULL, grid_by = 0.1,
                                    mc_n = 200) {
    n <- length(b0)
    if (is.null(tau)) tau <- rep(NA_real_, n)
    mc_n <- as.integer(mc_n)

    # Amplitude on link scale
    amplitude_link <- sqrt(b_sin^2 + b_cos^2)

    # Phase is reported as the fitted peak time in hours.
    phase_rad <- atan2(b_sin, b_cos)
    phase <- (phase_rad / (2 * pi) * period) %% period

    # MESOR and amplitude on response scale via numerical grid. In mixed mode,
    # average over U ~ N(0, tau^2) using deterministic normal quantile nodes.
    tgrid <- seq(0, period, by = grid_by)
    tgrid <- tgrid[tgrid < period]
    omega <- 2 * pi / period
    S <- sin(omega * tgrid)
    C <- cos(omega * tgrid)
    z <- stats::qnorm((seq_len(mc_n) - 0.5) / mc_n)

    amplitude_response <- rep(NA_real_, n)
    mesor_response <- rep(NA_real_, n)

    for (i in seq_len(n)) {
        if (is.na(b0[i]) || is.na(b_sin[i]) || is.na(b_cos[i])) next
        eta <- b0[i] + b_sin[i] * S + b_cos[i] * C
        if (is.finite(tau[i]) && tau[i] > 0) {
            mu <- colMeans(.invlogit(outer(tau[i] * z, eta, "+")))
        } else {
            mu <- .invlogit(eta)
        }
        mu_bar <- mean(mu)
        a_mu <- (max(mu) - min(mu)) / 2
        amplitude_response[i] <- a_mu
        mesor_response[i] <- mu_bar
    }

    list(
        amplitude_link      = amplitude_link,
        amplitude_response  = amplitude_response,
        phase               = phase,
        mesor_link          = b0,
        mesor_response      = mesor_response
    )
}

#' Return top CpGs by significance
#'
#' Filters to converged CpGs, selects key columns, and returns the top
#' `n` sorted by `sort_by`.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] object
#'   returned by [fitCosinor()].
#' @param n Integer, number of top CpGs to return (default 10).
#' @param sort_by Character, column to sort by (default `"pvalue"`).
#' @return A [S4Vectors::DataFrame] with the top CpGs.
#' @export
#'
#' @examples
#' ex <- readRDS(system.file("extdata", "example_data.rds",
#'                          package = "RhyME"))
#' se <- ex$se[seq_len(5), ]
#' fit <- fitCosinor(se, time_col = "time", subject_col = "subject",
#'                   prior = "none", verbose = FALSE)
#' topCpGs(fit, n = 3)
topCpGs <- function(se, n = 10, sort_by = "pvalue") {
    res <- rowData(se)
    conv <- res$converged
    conv[is.na(conv)] <- FALSE

    sub <- res[conv, ]

    out <- DataFrame(
        cpg                = rownames(sub),
        amplitude_response = sub$amplitude_response,
        phase              = sub$phase,
        mesor_response     = sub$mesor_response,
        pvalue             = sub$pvalue,
        padj               = sub$padj
    )

    ord <- order(out[[sort_by]], na.last = TRUE)
    out <- out[ord, ]
    n <- min(n, nrow(out))
    out[seq_len(n), ]
}
