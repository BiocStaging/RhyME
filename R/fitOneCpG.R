#' Fit a single CpG with the compiled TMB beta-regression cosinor model
#'
#' Fits the full (sin + cos) and null (intercept) models with the package's compiled TMB
#' template (`beta_glmm` for repeated-measures designs, `beta_fe` for
#' fixed / cross-sectional designs) and derives either a joint Wald test or a
#' likelihood-ratio test. Replaces the former glmmTMB backend.
#'
#' @param cpg_id Character, the CpG identifier.
#' @param y Numeric vector of beta values for this CpG.
#' @param df A data.frame with columns: `sin_term`, `cos_term`, and optionally
#'   `subject`. Must be aligned with `y`.
#' @param X_cov Numeric design matrix for optional adjustment covariates.
#' @param mode Character, either `"fixed"` or `"mixed"`.
#' @param prior_phi Named list with `mean` and `sd` for the normal prior on
#'   log(phi), or NULL for no prior (pilot pass / EB disabled).
#' @param eps Numeric, squeeze bound for beta values.
#' @param test Character, one of `"none"`, `"Wald"`, or `"LRT"`.
#'
#' @return A named list containing coefficients, dispersion, the selected
#'   test result, and convergence status.
#' @noRd
.fitOneCpG <- function(cpg_id, y, df, X_cov = NULL, mode = "fixed",
                       prior_phi = NULL, eps = 1e-4,
                       test = c("none", "Wald", "LRT")) {

    test <- match.arg(test)

    na_result <- list(
        cpg = cpg_id,
        b0 = NA_real_, b_sin = NA_real_, b_cos = NA_real_,
        beta = numeric(),
        log_phi = NA_real_, tau = NA_real_,
        pvalue = NA_real_,
        converged = FALSE
    )

    # Prepare data
    df$beta_val <- y
    keep <- !is.na(df$beta_val) & is.finite(df$beta_val)
    df <- df[keep, , drop = FALSE]
    if (is.null(X_cov)) {
        X_cov <- matrix(numeric(), nrow = length(y), ncol = 0L)
    }
    X_cov <- as.matrix(X_cov)[keep, , drop = FALSE]

    df$beta_val <- .squeezeValues(df$beta_val, c(eps, 1 - eps))
    yv <- df$beta_val

    # Dispersion prior (EB shrinkage). use_prior = 0 -> plain MLE.
    use_prior <- if (is.null(prior_phi)) 0L else 1L
    pm <- if (use_prior) prior_phi$mean else 0
    ps <- if (use_prior) prior_phi$sd   else 1

    # Covariates are present in both models; rhythmicity adds sin and cos.
    design <- .makeDesignMatrices(df, X_cov)
    X1 <- design$full
    X0 <- design$null
    circ_idx <- match(c("sin_term", "cos_term"), colnames(X1))

    # Moment-based starting value for log(phi)
    mu0 <- mean(yv)
    v0 <- var(yv)
    phi0 <- mu0 * (1 - mu0) / v0 - 1
    lphi_init <- if (use_prior) pm else log(max(phi0, 1))

    if (mode == "mixed") {
        id0 <- as.integer(factor(df$subject)) - 1L
        nid <- max(id0) + 1L
    }

    .fitTMB <- function(X) {
        if (mode == "mixed") {
            dat <- list(model = "beta_glmm", y = yv, X = X, id = id0,
                        prior_mean = pm, prior_sd = ps, use_prior = use_prior)
            par <- list(beta = c(stats::qlogis(mu0), rep(0, ncol(X) - 1)),
                        log_phi = lphi_init, log_tau = log(0.1),
                        u = rep(0, nid))
            obj <- TMB::MakeADFun(dat, par, random = "u",
                                  DLL = "RhyME", silent = TRUE)
        } else {
            dat <- list(model = "beta_fe", y = yv, X = X,
                        prior_mean = pm, prior_sd = ps, use_prior = use_prior)
            par <- list(beta = c(stats::qlogis(mu0), rep(0, ncol(X) - 1)),
                        log_phi = lphi_init)
            obj <- TMB::MakeADFun(dat, par, DLL = "RhyME", silent = TRUE)
        }
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
        list(obj = obj, opt = opt, ll = -opt$objective)
    }

    tryCatch({
        m1 <- .fitTMB(X1)
        m0 <- NULL
        pvalue <- NA_real_

        # Likelihood-ratio test (full vs null), df = 2.
        if (test == "LRT") {
            m0 <- .fitTMB(X0)
            statistic <- 2 * (m1$ll - m0$ll)
            pvalue <- stats::pchisq(
                statistic, df = 2, lower.tail = FALSE
            )
        }

        # Fixed-effect coefficients; the final two columns are sin and cos.
        bidx <- which(names(m1$opt$par) == "beta")
        bhat <- stats::setNames(m1$opt$par[bidx], colnames(X1))

        # Covariance-aware joint Wald test for H0: b_sin = b_cos = 0.
        # The inverse Hessian includes covariance with the nuisance parameters;
        # its beta block is then used without assuming that the sine and cosine
        # coefficient estimates are independent.
        if (test == "Wald") {
            pvalue <- tryCatch({
                H <- stats::optimHess(m1$opt$par, m1$obj$fn, m1$obj$gr)
                H <- (H + t(H)) / 2
                vc <- solve(H)
                vc <- (vc + t(vc)) / 2
                beta_vcov <- vc[bidx, bidx, drop = FALSE]

                circ_vcov <- beta_vcov[circ_idx, circ_idx, drop = FALSE]
                circ_beta <- unname(bhat[circ_idx])
                circ_vcov_inv <- chol2inv(chol(circ_vcov))
                stat <- drop(crossprod(circ_beta,
                                       circ_vcov_inv %*% circ_beta))
                if (!is.finite(stat)) {
                    stop("invalid joint Wald statistic")
                }
                stat <- max(stat, 0)
                stats::pchisq(stat, df = 2, lower.tail = FALSE)
            }, error = function(e) NA_real_)
        }

        # Dispersion and random-intercept SD
        log_phi <- unname(m1$opt$par[names(m1$opt$par) == "log_phi"])
        tau <- NA_real_
        if (mode == "mixed") {
            tau <- exp(unname(m1$opt$par[names(m1$opt$par) == "log_tau"]))
        }

        converged <- m1$opt$convergence == 0 &&
            (test != "LRT" || m0$opt$convergence == 0)

        list(
            cpg = cpg_id,
            b0 = unname(bhat[1]),
            b_sin = unname(bhat["sin_term"]),
            b_cos = unname(bhat["cos_term"]),
            beta = bhat,
            log_phi = log_phi,
            tau = tau,
            pvalue = pvalue,
            converged = converged
        )
    }, error = function(e) {
        na_result
    })
}
