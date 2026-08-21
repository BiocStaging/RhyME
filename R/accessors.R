#' Compute circadian parameters from cosinor coefficients
#'
#' @param b0 Numeric vector of intercepts.
#' @param b_sin Numeric vector of sin coefficients.
#' @param b_cos Numeric vector of cos coefficients.
#' @param period Numeric, the period in hours.
#' @param tau Numeric vector of random-intercept SD estimates. Used to
#'   recover marginal response-scale curves in mixed mode.
#' @return A list with amplitude_link, amplitude_response, phase,
#'   mesor_link, mesor_response.
#' @noRd
.computeCircadianParams <- function(b0, b_sin, b_cos, period = 24,
                                    tau = NULL) {
    n <- length(b0)
    if (is.null(tau)) tau <- rep(NA_real_, n)

    # Amplitude on link scale
    amplitude_link <- sqrt(b_sin^2 + b_cos^2)

    # Phase is reported as the fitted peak time in hours.
    phase_rad <- atan2(b_sin, b_cos)
    phase <- (phase_rad / (2 * pi) * period) %% period

    amplitude_response <- rep(NA_real_, n)
    mesor_response <- rep(NA_real_, n)

    valid <- is.finite(b0) & is.finite(b_sin) & is.finite(b_cos)
    valid_idx <- which(valid)
    if (length(valid_idx)) {
        tau_use <- ifelse(is.finite(tau) & tau > 0, tau, 0)
        gh <- .normalGH(15L)

        # The marginal inverse-logit is monotone, so the response-scale extrema
        # occur exactly at b0 +/- amplitude_link; no time grid is needed here.
        mu_high <- numeric(length(valid_idx))
        mu_low <- numeric(length(valid_idx))
        for (k in seq_along(gh$nodes)) {
            shift <- tau_use[valid_idx] * gh$nodes[k]
            mu_high <- mu_high + gh$weights[k] * stats::plogis(
                b0[valid_idx] + amplitude_link[valid_idx] + shift
            )
            mu_low <- mu_low + gh$weights[k] * stats::plogis(
                b0[valid_idx] - amplitude_link[valid_idx] + shift
            )
        }
        amplitude_response[valid_idx] <- (mu_high - mu_low) / 2

        # Average the marginal curve over one period using 48 equally spaced
        # phase nodes and block vectorization across CpGs.
        theta <- 2 * pi * (0:47) / 48
        sin_t <- sin(theta)
        cos_t <- cos(theta)
        starts <- seq.int(1L, length(valid_idx), by = 1000L)
        for (start in starts) {
            idx <- valid_idx[
                start:min(start + 999L, length(valid_idx))
            ]
            eta <- outer(b0[idx], rep(1, 48L)) +
                outer(b_sin[idx], sin_t) + outer(b_cos[idx], cos_t)
            marginal_mu <- matrix(0, nrow = length(idx), ncol = 48L)
            for (k in seq_along(gh$nodes)) {
                marginal_mu <- marginal_mu + gh$weights[k] *
                    stats::plogis(eta + tau_use[idx] * gh$nodes[k])
            }
            mesor_response[idx] <- rowMeans(marginal_mu)
        }
    }

    list(
        amplitude_link      = amplitude_link,
        amplitude_response  = amplitude_response,
        phase               = phase,
        mesor_link          = b0,
        mesor_response      = mesor_response
    )
}

#' Gauss-Hermite nodes and weights for expectations under N(0, 1)
#' @noRd
.normalGH <- function(n) {
    J <- matrix(0, nrow = n, ncol = n)
    idx <- seq_len(n - 1L)
    off_diag <- sqrt(idx / 2)
    J[cbind(idx, idx + 1L)] <- off_diag
    J[cbind(idx + 1L, idx)] <- off_diag
    eig <- eigen(J, symmetric = TRUE)
    ord <- order(eig$values)
    list(
        nodes = sqrt(2) * eig$values[ord],
        weights = eig$vectors[1L, ord]^2
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
