# Rhythm summaries, marginal over the fitted random intercept.
.computeCircadianParams <- function(b0, b_sin, b_cos, period = 24,
                                    tau = NULL) {
    n <- length(b0)
    if (is.null(tau)) tau <- rep(NA_real_, n)

    amplitude_link <- sqrt(b_sin^2 + b_cos^2)

    # Phase is reported as the fitted peak time in hours.
    phase_rad <- atan2(b_sin, b_cos)
    phase <- (phase_rad / (2 * pi) * period) %% period
    phase[!is.finite(amplitude_link) | amplitude_link == 0] <- NA_real_

    amplitude_response <- rep(NA_real_, n)
    mesor_response <- rep(NA_real_, n)

    valid <- is.finite(b0) & is.finite(b_sin) & is.finite(b_cos) &
        is.finite(tau) & tau >= 0
    valid_idx <- which(valid)
    if (length(valid_idx)) {
        gh <- .normalGH(15L)

        # Monotonicity places the marginal extrema at b0 +/- amplitude_link.
        mu_high <- numeric(length(valid_idx))
        mu_low <- numeric(length(valid_idx))
        for (k in seq_along(gh$nodes)) {
            shift <- tau[valid_idx] * gh$nodes[k]
            mu_high <- mu_high + gh$weights[k] * stats::plogis(
                b0[valid_idx] + amplitude_link[valid_idx] + shift
            )
            mu_low <- mu_low + gh$weights[k] * stats::plogis(
                b0[valid_idx] - amplitude_link[valid_idx] + shift
            )
        }
        amplitude_response[valid_idx] <- (mu_high - mu_low) / 2

        # Average over one period, processing CpGs in blocks to limit memory.
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
                    stats::plogis(eta + tau[idx] * gh$nodes[k])
            }
            mesor_response[idx] <- rowMeans(marginal_mu)
        }
    }

    list(
        logit_amplitude = amplitude_link,
        amplitude = amplitude_response,
        phase = phase,
        mesor = mesor_response
    )
}

# Gauss-Hermite quadrature for N(0, 1).
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
